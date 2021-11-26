pacman::p_load("devtools", "RcppArmadillo", "Rcpp", "reticulate", "rjson", "BNPmix")
library("devtools")
library(RcppArmadillo)
library(Rcpp)
library("reticulate")
library(rjson)
setwd("~/Documents/git/msMK")
source("tests/tests.R")
devtools::load_all()

setwd("~/Documents/git/msMK/tests/estimation/")
opts = fromJSON(file="settings.json")
mod_sim = opts$mod_sim

# library(DPpackage)
library(BNPmix)

BNPmix.cpo = function(mcmc){
  cpo = colMeans(1/mcmc$density)^(-1)
  return(cpo)
}

BNPmix.lpml = function(mcmc){
  mean(log(BNPmix.cpo(mcmc)))
}

lpml_mat = matrix(nrow = mod_sim, ncol = 2)
for(i in 1:mod_sim){
  delta = opts$delta
  indep = opts$indep
  nsim = opts$nsim
  burnin = opts$burnin
  p = as.integer(opts$p)
  n = as.integer(opts$n)

  if(p > 2){
    gen.str = paste0(opts$gen, "_p" )
    df.str = paste0(opts$df, "_p")
    gen = function(n) get(gen.str)(n, p = p)
    df = function(x) get(df.str)(x, p = p)
  } else if (p == 2){
    gen.str = paste0(opts$gen, "_2" )
    df.str = paste0(opts$df, "_2")
    gen = get(gen.str)
    df = get(df.str)
  }
  
  dir.str = paste0("~/Documents/git/msMK/tests/estimation/", opts$gen, "_", p)
  if(!dir.exists(dir.str)){
    dir.create(dir.str)
    dir.create(paste0(dir.str, "/img"))
    dir.create(paste0(dir.str, "/tex"))
  }
  setwd(dir.str)

  cat("Simulazione -- ", i, "\n")
  set.seed(i)
  
  y = gen(n) 
  y.stand = y
  
  
  name = paste0("i", i,"-PYfix")

  # Initial state
  state <- NULL

  # MCMC parameters
  mu0 = apply(y, 2, mean)
  sig0 = diag(apply(y, 2, var))

  nburn <- burnin
  nsave <- nsim + burnin
  nskip <- 0
  ndisplay <- 100
  mcmc <- list(nburn=nburn,niter=nsave,method="ICS",hyper=FALSE)

  # PY without hyperprior
  prior <- list(m0 = mu0, Sigma0 = sig0, n0=p+2)

  grid <- expand.grid(seq(-7, 7, length.out = 50),
                      seq(-7, 7, length.out = 50))
  # PY plus hyperprior
  # prior <- list(a0=1,b0=1, m2=mu0, psiinv2=sig0,
  #                nu1=p, nu2=p, 
  #                s2 = diag(1, nrow=p),
  #                tau1=2,tau2=200)


  fit1.1 <- NULL
  attempt <- 1
  while( is.null(fit1.1) && attempt <= 50 ) {
    attempt <- attempt + 1
    try(
        fit1.1 <- PYdensity(y, prior=prior,mcmc=mcmc, output=list(grid=data.frame(y)))
    )
  } 

    
  cat("Done", "\n\n")
  
  cat("Calculating LPML", "\n")
  lpml_dp = BNPmix.lpml(fit1.1)
  lpml_true = lpml(y.stand, df)
  lpml_mat[i, ] = c(lpml_dp, lpml_true)
  colnames(lpml_mat) = c("PY", "True")
  cat("Done", "\n\n")
}
save(lpml_mat, file=paste0("lpml-PYfix-n", n,".Rdata"))
lpml_mean = as.matrix(apply(lpml_mat, 2, mean, na.rm = TRUE))
lpml_sd = as.matrix(apply(lpml_mat, 2, sd, na.rm = TRUE))
library(magrittr)
library(kableExtra)
tab_name = paste0("lpml-PYfix-n", n)
cbind(lpml_mean, lpml_sd) %>%
  set_colnames(c("lpml", "sd")) %>%
  kbl("latex", booktabs = TRUE, escape = FALSE,
      col.names = c("$\\overbar{\\text{lpml}}$", "$\\text{sd}(\\overbar{\\text{lpml}})$"),
      caption = "Average lpml and related standard error over $n_{\\text{sim}} = 10$ simulations.",
      label = ) %>%
  kable_styling(latex_options = "HOLD_position") %>%
  cat(file=paste0("tex/", tab_name, ".tex"))

# # Initial state
# state <- NULL

# # MCMC parameters

# nburn <- 1000
# nsave <- 400
# nskip <- 0
# ndisplay <- 100
# mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)

# # Example of Prior information 1
# # Fixing alpha, m1, and Psi1

# prior1 <- list(alpha=1,m1=mu0,psiinv1=sig0,nu1=4,
#                tau1=1,tau2=4)

# fit1.1 <- DPdensity(y=y, ngrid = 2500, prior=prior1,mcmc=mcmc, state=state,status=TRUE)

# contour(fit1.1$x1, fit1.1$x2, fit1.1$dens, ylim = c(-4, 4), xlim = c(-4,4))
