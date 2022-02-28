#' Read thresholds from bounds.csv file and return it as a
#' list of lower and upper bounds
readThresholds = function(){
  tbl <- as.matrix(read.csv("bounds.csv", header = FALSE))
  n = NROW(tbl)/2
  lb <- matrix(NA, ncol = p, nrow = n)
  ub <- matrix(NA, ncol = p, nrow = n)
  for(i in 1:n){
    lb[i, ] <- tbl[2*i - 1 , ]
    ub[i, ] <- tbl[2*i, ]
  }
  return(list("lb" = lb, "ub" = ub))
}

#' Sample from the msMK mixture density
#'
#' @param n: number of observations to generate 
#' @param TH: location parameters for the density 
#' @param SIG: sequence of scale parameters for the density
#' @param prob: probabilities of each cluster
#'
#' @return a matrix of size n x p containing the generated data 
rmsMK = function(n, mcmc, idx){
  require(mvtnorm)
  TH = mcmc$theta[[idx]]
  SIG = mcmc$sigma[[idx]]
  prob = mcmc$prob[[idx]]
  out = matrix(0, nrow = n, ncol = length(TH[[1]]))
  clus = sample(1:length(prob), size = n, prob = prob, replace = TRUE)
  for(i in 1:length(TH)){
    clus_i = (clus == i)
    n_i = sum(clus_i)
    if(n_i != 0){
      theta = as.vector(TH[[i]])
      sigma = SIG[[i]]
      out[clus_i, ] = mvtnorm::rmvnorm(n_i, theta, sigma)
    }
  }
  return(out)
}

#' Sample from the msMK mixture density
#'
#' @param n: number of observations to generate 
#' @param TH: location parameters for the density 
#' @param SIG: sequence of scale parameters for the density
#' @param prob: probabilities of each cluster
#'
#' @return a matrix of size n x p containing the generated data 
posterior_predict_msMK = function(mcmc, y, ndraws = NULL){
  require(mvtnorm)
  n = NROW(y)
  p = NCOL(y)
  nsim = length(mcmc$theta)
  burnin = mcmc$burnin
  thrs  = mcmc$thrs
  if(is.null(ndraws)){
    # Draw the whole dataset for each sampled posterior parameter
    out = array(0, dim = c(n, p, nsim - burnin))
    for(b in (burnin+1):nsim){
      TH = mcmc$theta[[b]]
      SIG = mcmc$sigma[[b]]
      prob = mcmc$prob[[b]]
      idx = sample(1:length(prob), size = n, prob = prob, replace = TRUE)
      
      for(i in 1:length(prob)){
        idx_i = (idx == i)
        n_i = sum(idx_i)
        if(n_i != 0){
          theta = as.vector(TH[, i])
          sigma = SIG[, , i]
          out[idx_i, , b - burnin] = mvtnorm::rmvnorm(n_i, theta, sigma)
        }
      }
    }
  } else{
    out = array(0, dim = c(n, p, ndraws))
    pts = sample((burnin+1):nsim, size = ndraws, replace = TRUE)
    for(j in 1:length(pts)){
      b = pts[j]
      TH = mcmc$theta[[b]]
      SIG = mcmc$sigma[[b]]
      prob = mcmc$prob[[b]]
      idx = sample(1:length(prob), size = n, prob = prob, replace = TRUE)
      
      for(i in 1:length(prob)){
        idx_i = (idx == i)
        n_i = sum(idx_i)
        if(n_i != 0){
          theta = as.vector(TH[, i])
          sigma = SIG[ , , i]
          out[idx_i, , j] = mvtnorm::rmvnorm(n_i, theta, sigma)
        }
      }
    }
  }
  return(out)
}


#' Mean posterior predictive distribution for a fitted msMK model
#'
#' @param n: number of observations to generate 
#' @param TH: location parameters for the density 
#' @param SIG: sequence of scale parameters for the density
#' @param prob: probabilities of each cluster
#'
#' @return a matrix of size n x p containing the generated data 
dmsMK_ave = function(x, mcmc){
  require(TruncatedNormal)
  out = rep(0, NROW(x))
  nsim = length(mcmc$theta)
  burnin = mcmc$burnin
  for(b in (burnin+1):nsim){
    out = out + dmsMK(x, mcmc$theta[[b]], mcmc$sigma[[b]], mcmc$thrs, mcmc$prob[[b]], mcmc$indep)
  }
  out = out / (nsim - burnin)
  return(out)
}


#' Calculate logarithm of the pseudo-marginal likelihood (LPML) for a set of observations x 
#' under the multiscale mixture of kernels
#' 
#' @param x: matrix of observed data
#' @param mcmc: list of realizations from the posterior distribution of the data (without burnin) 
#'
#' @return a single value, the estimated LPML 
msMK.cpo = function(x, mcmc){
  require(TruncatedNormal)
  B = length(mcmc$theta)
  burnin = mcmc$burnin
  CPO_inv = rep(0, length = NROW(x))       # Create vector of CPO_i ^{-1}
  thrs  = mcmc$thrs
  indep = mcmc$indep
  for(b in (burnin+1):B){
    print(b)
    theta = mcmc$theta[[b]]
    sigma = mcmc$sigma[[b]]
    prob = mcmc$prob[[b]]
    CPO_inv = CPO_inv + (1 / as.vector(dmsMK(x, theta, sigma, thrs, prob, indep)))
  }
  CPO = (CPO_inv / (B - burnin))^(-1)
  return(CPO)
}

msMK.lpml = function(cpo){
  return(mean(log(cpo)))
}

#' Calculate logarithm of the pseudo-marginal likelihood (LPML) for a set of observations x 
#' under a known density function
#' 
#' @param x: matrix of observed data
#' @param df: density function
#' @param ...: further arguments to df
#'
#' @return a single value, the estimated LPML 
lpml = function(x, df, ...){
  CPO = NULL
  attempt <- 50
  while( is.null(CPO) && attempt <= 50 ) {
    attempt <- attempt + 1
    try(
        CPO <- df(x, ...)
    )
  } 
  return(mean(log(CPO)))
}


#' Contour plot of the predictive density for a single parameter draw from the
#' posterior msMK distribution
#' 
#' @param mcmc: output of the msMK_mcmc function 
#' @param idx: integer, index of the parameter from which to calculate the density
#' @param npoints: number of points in each dimension
#' @param y: provide the data matrix in order to superimpose the observed points to the contour plot
#' @param ...: futher arguments to contour or filled.contour
#'
#' @return a list containing the sequence of points in each dimension, along with
#' the calculated density values 
plot_single_2d = function(mcmc, idx, npoints = 50, y = NULL, ...){
  theta = mcmc$theta[[idx]]
  sigma = mcmc$sigma[[idx]]
  prob = mcmc$prob[[idx]]
  thrs = mcmc$thrs
  
  if(!is.null(y)){
    xx = seq(min(y[ , 1]), max(y[ , 1]), length.out = npoints)
    yy = seq(min(y[ , 2]), max(y[ , 2]), length.out = npoints)
    grid = as.matrix(expand.grid(xx,yy))
    zz = matrix(dmsMK(grid, theta, sigma, thrs, prob), nrow = npoints)
    
    plot(y, pch = 20, col = "lightgray")
    contour(xx, yy, zz, add=TRUE, ...)
  } else{
    xx = seq(-4, 4, length.out = npoints)
    yy = seq(-4, 4, length.out = npoints)
    grid = as.matrix(expand.grid(xx,yy))
    zz = matrix(dmsMK(grid, theta, sigma, thrs, prob, mcmc$indep), nrow = npoints)
    
    filled.contour(xx, yy, zz, ...)
  }
  return(list("dim1" = xx, "dim2" = yy, "z" = zz))
}

#' Contour plot of the average predictive density from the posterior msMK distribution
#'
#' @param mcmc: output of the msMK_mcmc function 
#' @param npoints: number of points in each dimension
#' @param y: provide the data matrix in order to superimpose the observed points to the contour plot
#' @param ...: futher arguments to contour or filled.contour
#'
#' @return a list containing the sequence of points in each dimension, along with
#' the calculated density values 
plot_average_2d = function(mcmc, npoints = 50, y, addy = FALSE, ...){
  xx = seq(min(y[ , 1]), max(y[ , 1]), length.out = npoints)
  yy = seq(min(y[ , 2]), max(y[ , 2]), length.out = npoints)
  grid = as.matrix(expand.grid(xx,yy))
  zz = matrix(dmsMK_ave(grid, mcmc), nrow = npoints)
  if(addy){
    plot(y, pch = 20, col = "lightgray")
    contour(xx, yy, zz, add=TRUE, ...)
  } else{
    contour(xx, yy, zz, ...)
  }
  return(list("dim1" = xx, "dim2" = yy, "z" = zz))
}


predictive_smooth_average = function(mcmc, y, ylim = NULL, ndraws = NULL){
  require(KernSmooth)
  pred = posterior_predict_msMK(mcmc, y, ndraws)
  
  # Bandwidth di default smoothScatter
  bandwidth <- diff(apply(pred, 2, stats::quantile,
                          probs = c(0.05, 0.95),
                          na.rm = TRUE, names = FALSE)) / 25
  bandwidth[bandwidth==0] <- 1
  if(!is.null(ylim)){
    predSmooth = apply(pred, 3, bkde2D, bandwidth = bandwidth,
                       range.x = list(c(ylim[1], ylim[2]), c(ylim[1], ylim[2])))
  } else{
    predSmooth = apply(pred, 3, bkde2D, bandwidth = bandwidth,
                       range.x = list(c(min(y), max(y)), c(min(y), max(y))))
  }
  out = matrix(0, nrow = NROW(predSmooth[[1]]$fhat), ncol = NCOL(predSmooth[[1]]$fhat))
  for(i in 1:length(predSmooth)){
    out = out + predSmooth[[i]]$fhat
  }
  out = out / length(predSmooth)
  return(list("x1" = predSmooth[[1]]$x1, "x2" = predSmooth[[1]]$x2, "fhat" = out))
}


#' Plot posterior binary tree of the multiscale model
#'
#' @param mcmc: output of the msMK_mcmc function 
#' @param labels: type of labels to provide, `ordered` uses 1:(2^(smax+1)-1), whereas `sh` uses the (s,h) notation.
#'
#' @return NULL 
msMK.plot = function(mcmc, labels = c("ordered", "none", "sh"), ...){
  require(igraph)
  lab = match.arg(labels, c("ordered", "none", "sh"))
  nelem = length(mcmc$prob[[1]])
  B = length(mcmc$prob)
  burnin = mcmc$burnin
  G = graph.tree(nelem, children=2)
  co <- layout.reingold.tilford(G, params=list(root=1))
  
  # Calculate median of simulated data
  prob_list = mcmc$prob[(burnin+1):B]
  prob = apply(do.call(rbind, prob_list), 2, median)
  
  # Normalize medians to unit interval
  prob = (prob - min(prob) + min(prob)/nelem)/(max(prob) - min(prob) + min(prob)/nelem)
  
  # Vertex sizes proportional to normalized values
  size = 11 * prob
  
  if(lab == "sh"){
    # Use (s, h), s = 0, ..., smax, h = 1, ..., 2^s labels
    labels = NULL
    smax = log2(nelem + 1) - 1
    for(s in 0:smax){
      labels = c(labels, paste0(s, ",", 1:2^s))
    }
  } else if (lab == "none"){
    labels = NA
  }
  else{
    labels = 1:nelem
  }
  
  # Plot binary tree with sensible defaults
  plot(G, layout=co, vertex.size = size,
       vertex.label = labels,
       edge.width = 1,
       edge.arrow.width = rep(0, length(size)),
       vertex.size2 = 0,
       vertex.label.cex = 1,
       asp = 0.45,
       margin = -0.1,
       arrow.mode = 0,
       ...)
}

#' Plot expected tree depth as a function of MCMC draws
#'
#' @param mcmc: output of the msMK_mcmc function 
#' @param include.burnin: logical, if FALSE (default) then burnin simulations are not included
#'
#' @return NULL 
msMK.plot.depth = function(mcmc, include.burnin = FALSE){
  nelem = length(mcmc$prob[[1]])
  smax = log2(nelem + 1) - 1
  B = length(mcmc$prob)
  burnin = mcmc$burnin

  nsim = B - burnin

  # Get list of probabilities
  if(include.burnin){
    prob_list = mcmc$prob
  } else{
    prob_list = mcmc$prob[(burnin+1):B]
  }
  prob = do.call(rbind, prob_list)
  
  out = rep(0, length = NROW(prob))
  for(s in 0:smax){
    for(h in 1:2^s){
      out = out + s * prob[ , 2^s - 1 + h ]
    }
  }

  if(include.burnin){
    plot(1:NROW(prob), out, type = "l", ylim = c(0, smax),
         xlab = "sim", ylab = "E[S|-]")
    abline(v = burnin, lty = "dashed")
  } else{
    plot((burnin+1):B, out, type = "l", ylim = c(0, smax),
         xlab = "sim", ylab = "E[S|-]")
  }
}
