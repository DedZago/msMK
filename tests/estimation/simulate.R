library("devtools")
# install_github("DedZago/msMK")
library("msMK")
setwd("~/Documents/git/msMK")
library(RcppArmadillo)
library(Rcpp)
library("reticulate")
library(rjson)
source("tests/tests.R")
load_all()

reticulate::py_install(c("numpy", "numpy-hilbert-curve"), pip = TRUE)
reticulate::source_python("inst/hilbertSplit.py")

setwd("~/Documents/git/msMK/tests/estimation/")
opts = fromJSON(file="settings.json")
mod_sim = opts$mod_sim

lpml_mat = matrix(nrow = mod_sim, ncol = 2)

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
  
for(i in 1:mod_sim){
  # Skip if already saved model
  name = paste0("i", i,"-ind", indep,"-n",n)
  if(file.exists(paste0(name, ".Rdata"))){
    print(i)
    load(paste0(name, ".Rdata"))
    next
  }
  
  a = deltaToAlpha(delta)[3]
  cat("delta:", delta, ", ", "alpha:", a, "\n")
  cat("Simulazione -- ", i, "\n")
  set.seed(i)
  
  y = gen(n) 
  y.stand = y
  
  # if(p == 2){
  #   png(paste0(dir.str, "/img/true.png"), width = 2200, height = 1200, res = 100) 
  #   densContour(df, y.stand)
  #   dev.off()
  # }
 
  # ---- Create splits
  smax = 6L
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
  
  
  # if(p == 2){
  #   cat("Drawing smoothed predictive density", "\n")
  #   png(paste0("img/", name, "_ppdSmooth.png"), width = 2200, height = 1200, res = 100) 
  #   pred_mean = predictive_smooth_average(mcmc.msmk, y = y.stand)
  #   contour(pred_mean$x1, pred_mean$x2, pred_mean$fhat)
  #   dev.off()
  # }
  
  cat("Done", "\n\n")
  
  cat("Calculating LPML", "\n")
  lpml_msmk = msMK.lpml(mcmc.msmk$cpo)
  lpml_true = lpml(y.stand, df)
  lpml_mat[i, ] = c(lpml_msmk, lpml_true)
  colnames(lpml_mat) = c("msMK", "True")
  cat("Done", "\n\n")
  save.image(file = paste0(name, ".Rdata"))
}
save(lpml_mat, file=paste0("lpml-ind", indep, "-n", n, ".Rdata"))
lpml_mean = as.matrix(apply(lpml_mat, 2, mean, na.rm = TRUE))
lpml_sd = as.matrix(apply(lpml_mat, 2, sd, na.rm = TRUE))
library(magrittr)
library(kableExtra)
tab_name = paste0("lpml-", "-ind", indep, "-n", n)
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
