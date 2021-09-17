library("devtools")
library(RcppArmadillo)
library(Rcpp)
library("reticulate")
setwd("~/Documents/git/MSc-thesis/code/msMK")
source("tests/tests.R")
devtools::load_all()

reticulate::py_install(c("numpy", "numpy-hilbert-curve"), pip = TRUE)
reticulate::source_python("hilbertSplit.py")

setwd("~/Documents/git/MSc-thesis/code/tests/estimation/r9_10/")

mod_sim = 10

delta = 0.25
a = deltaToAlpha(delta)[2]
lpml_mat = matrix(nrow = mod_sim, ncol = 2)
for(i in 1:mod_sim){
  cat("delta:", delta, ", ", "alpha:", a, "\n")
  cat("Simulazione -- ", i, "\n")
  set.seed(i)
  n = 1000
  p = 10L
  gen = function(n) r9_p(n, p)
  df = f9_p

  y <- NULL
  attempt_gen <- 1
  while( is.null(y) && attempt_gen <= 50 ) {
    attempt_gen <- attempt_gen + 1
    try(
        y <- gen(n) 
    )
  } 
  y.stand = y

  # ---- Create splits
  smax = 5L
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
  args = commandArgs(trailingOnly=TRUE)
  if(length(args) >= 1){
    if(args[1] == "FALSE"){
      indep = FALSE
    } else{
      indep = TRUE
    }
    print(indep)
  } else{
    stop("Independence not provided (TRUE/FALSE)")
  }

  nsim = 200
  burnin = 800

  mcmc.msmk <- NULL
  attempt <- 1
  while( is.null(mcmc.msmk) && attempt <= 50 ) {
    attempt <- attempt + 1
    try(
        mcmc.msmk <- msMK_mcmc_test(nsim, y.stand, a, b, delta, smax, mu0, k0, sig0, indep, thrs$lb, thrs$ub,
                                    burnin = burnin)
    )
  } 

  name = paste0("r9_10-d", delta, "-a", a, "-i", i, "-ind",indep)
  save.image(file = paste0(name, ".Rdata"))

  cat("Calculating LPML", "\n")
  lpml_msmk = msMK.lpml(y.stand, mcmc.msmk)
  lpml_true = lpml(y.stand, df)
  lpml_mat[i, ] = c(lpml_msmk, lpml_true)
  colnames(lpml_mat) = c("msMK", "True")
  cat("Done", "\n\n")
  save.image(file = paste0(name, ".Rdata"))
}
lpml_mean = as.matrix(apply(lpml_mat, 2, mean, na.rm = TRUE))
lpml_sd = as.matrix(apply(lpml_mat, 2, sd, na.rm = TRUE))
library(magrittr)
library(kableExtra)
tab_name = paste0("lpml-", "r9_10-d", delta, "-a", a, "-ind",indep)
cbind(lpml_mean, lpml_sd) %>%
  set_colnames(c("lpml", "sd")) %>%
  kbl("latex", booktabs = TRUE, escape = FALSE,
      col.names = c("$\\overbar{\\text{lpml}}$", "$\\text{sd}(\\overbar{\\text{lpml}})$"),
      caption = "Average lpml and related standard error over $n_{\\text{sim}} = 30$ simulations.",
      label = ) %>%
  kable_styling(latex_options = "HOLD_position") %>%
  cat(file=paste0("tex/", tab_name, ".tex"))
