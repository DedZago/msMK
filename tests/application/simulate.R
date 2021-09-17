library("devtools")
library(RcppArmadillo)
library(Rcpp)
library("reticulate")
library(rjson)
setwd("~/Documents/git/MSc-thesis/code/msMK")
source("tests/tests.R")
devtools::load_all()

reticulate::py_install(c("numpy", "numpy-hilbert-curve"), pip = TRUE)
reticulate::source_python("hilbertSplit.py")

setwd("~/Documents/git/MSc-thesis/code/tests/application/")

delta = 0.25
alpha = deltaToAlpha(delta)
a = alpha[3]
indep = FALSE
nsim = 2000
burnin = 1000
dati = read.csv("new_galaxy_data.csv")

colnames(dati) = c("colour", "luminosity", "log(density)", "redshift", "velocity.of.cluster")
y = dati[ , -5 ]
y.stand = y
y.stand = as.matrix(y.stand)

png("marginal.png", width = 1500, height = 900, res = 100) 
par(mfrow = c(2,2))
for(j in 1:4){
  hist(y[, j], nclass = 100, prob = TRUE, main="", xlab = colnames(dati)[j])
}
par(mfrow = c(1,1))
dev.off()

png("pairs.png", width = 1500, height = 900, res = 100) 
pairs(y)
  
  # if(p == 2){
  #   png(paste0(dir.str, "/img/true.png"), width = 2200, height = 1200, res = 100) 
  #   densContour(df, y.stand)
  #   dev.off()
  # }
 
  # ---- Create splits
p = NCOL(y.stand)
smax = 8L
hilbertSplit(p, smax)
thresholds = readThresholds()
b = 1
mu0 = apply(y, 2, mean) 
k0 = apply(y, 2, var)

thrs = list()
thrs$lb = apply(thresholds$lb, 2, qnorm, mean = mu0, sd = sqrt(k0))     # lower quantiles
thrs$ub = apply(thresholds$ub, 2, qnorm, mean = mu0, sd = sqrt(k0))     # upper quantiles

# sig0 = diag(p)
sig0 = diag(apply(y, 2, var))

mcmc.msmk <- NULL
attempt <- 1
while( is.null(mcmc.msmk) && attempt <= 50 ) {
attempt <- attempt + 1
  try(
      mcmc.msmk <- msMK_mcmc_test(nsim, y.stand, a, b, delta, smax, mu0, k0, sig0, indep, thrs$lb, thrs$ub,
                                  burnin = burnin)
  )
} 

cat("Saving model image\n\n")
save.image(paste0("model-estimated-", indep, ".Rdata"))
# load("model-estimated-FALSE.Rdata")

cat("Plot expected depth\n\n")
png(paste0("expected-depth-", indep, ".png"), width = 1500, height = 900, res = 100) 
msMK.plot.depth(mcmc.msmk, TRUE)
dev.off()

cat("Plot binary tree\n\n")
png(paste0("binary-tree-", indep, ".png"), width = 1500, height = 900, res = 100) 
msMK.plot(mcmc.msmk, "none")
dev.off()

set.seed(123)
npred = 500
pred = posterior_predict_msMK(mcmc.msmk, y.stand, ndraws = npred)

library(stringr)

temp = matrix(NA, nrow = NROW(mcmc.msmk$cluster), ncol = dim(mcmc.msmk$cluster)[3])
for(k in 1001:dim(mcmc.msmk$cluster)[3]){
  print(k)
  slice = mcmc.msmk$cluster[ , , k ]
  temp[, k] = apply(slice, 1, function(row) paste0(as.character(row[1]), as.character(row[2])))
}

cluster_cnt = apply(temp, 1, table)
cluster = unlist(lapply(cluster_cnt, function(l) names(l)[which.max(l)]))
cluster = factor(cluster)
num_colors <- nlevels(cluster)
color <- colorRampPalette(num_colors)

predictive_fig10_plot = function(mcmc, y, pred, cols, nlevels = 20){
  require(KernSmooth)

  colnames(y) = c("color","luminosity","log(density)","redshift")
  # Bandwidth di default smoothScatter
  if(cols[1] == cols[2]){
    npred = dim(pred)[3]
    j = cols[1]
    hist(dati[, j], nclass = 100, prob = TRUE, main=colnames(y)[j], xlab = colnames(dati)[j])
    xlim = c(min(dati[,j]), max(dati[,j]))
    for(i in 1:npred){
      temp = pred[ , j, i]
      temp = temp[temp < max(y.stand[,j]) + 6*sd(y.stand[, j])]
      temp = temp[temp > min(y.stand[,j]) - 6*sd(y.stand[, j])]
      lines(density(temp, from = xlim[1], to = xlim[2]), col = "blue")
    }

  } else{
    pred = pred[ , cols, ]
    y = y[, cols]
    bw <- diff(apply(pred, 2, stats::quantile,
                            probs = c(0.05, 0.95),
                            na.rm = TRUE, names = FALSE)) / 25
    bw[bw==0] <- 1
    lims = list(c(min(y[ , 1 ]), max(y[, 1])),
                c(min(y[ , 2 ]), max(y[, 2])))
    predSmooth = apply(pred, 3, bkde2D, bandwidth = bw,
                       range.x = lims)
    
    out = matrix(0, nrow = NROW(predSmooth[[1]]$fhat), ncol = NCOL(predSmooth[[1]]$fhat))
    for(i in 1:length(predSmooth)){
      out = out + predSmooth[[i]]$fhat
    }
    out = out / length(predSmooth)
    plot(y, pch = 20, cex = 0.1, xlab = colnames(y)[1], ylab = colnames(y)[2], col = cluster)
    contour(predSmooth[[1]]$x1, predSmooth[[1]]$x2, out,
            nlevels = nlevels, add = TRUE, col = "black")
    return(list("x1" = predSmooth[[1]]$x1, "x2" = predSmooth[[1]]$x2, "fhat" = out, "y" = y))
  }
}

png(paste0("posterior-predict-", indep, ".png"), width = 1500, height = 900, res = 100) 
par(mfrow = c(2,3))
for(i in 1:3){
  for(k in (i+1):4){
    if(i == 3 | k == 3){
      nlevels = 15
    } else{
      nlevels = 15
    }
    predictive_fig10_plot(mcmc.msmk, y, pred, c(k,i), nlevels)
  }
}
par(mfrow = c(1,1))
dev.off()

# Calculate LPML from model
cpo = msMK.cpo(y.stand, mcmc.msmk)

save(cpo, file=paste0("cpo-", indep, ".Rdata"))
lpml = msMK.lpml(cpo)
print(lpml)
