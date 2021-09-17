library("devtools")
library(RcppArmadillo)
library(Rcpp)
library("reticulate")
setwd("~/Documents/git/MSc-thesis/code/msMK")
source("tests/tests.R")
devtools::load_all(recompile = TRUE)

reticulate::py_install(c("numpy", "numpy-hilbert-curve"), pip = TRUE)
reticulate::source_python("hilbertSplit.py")

setwd("~/Documents/git/MSc-thesis/code/tests/estimation")

# ---- Generate data
set.seed(123)
n = 10000
p = 2L
gen = r5_2
df = f5_2

y = gen(n, p) 
# y.stand = scale(y)
y.stand = y

# ---- Create splits
smax = 5L
hilbertSplit(p, smax)
thresholds = readThresholds()
a = 1
b = 1
delta = 0.1
mu0 = rep(0, p)
# k0 = rep(1, p)
k0 = apply(y, 2, var)

thrs = list()
thrs$lb = apply(thresholds$lb, 2, qnorm, sd = sqrt(k0))     # lower quantiles
thrs$ub = apply(thresholds$ub, 2, qnorm, sd = sqrt(k0))     # upper quantiles

sig0 = diag(p)
indep = FALSE

nsim = 200
burnin = 300
mcmc.msmk = msMK_mcmc_test(nsim, y.stand, a, b, delta, smax, mu0, k0, sig0, indep, thrs$lb, thrs$ub,
                           burnin = burnin)

save.image(file = paste0("r5_2-a", a, "-b", b,"-d", delta, "-S", smax,".Rdata"))

#---- Contour for a single parameter update
smoothScatter(y.stand)
plot_single_2d(mcmc.msmk, nsim + burnin)
plot_single_2d(mcmc.msmk, nsim + burnin, y = y.stand)

smoothScatter(y.stand)
plot_single_2d(mcmc.msmk, nsim + burnin)
plot_average_2d(mcmc.msmk, y = y.stand)
{
  par(mfrow = c(1,2))
  smoothScatter(y.stand, ylim = c(-3, 3), xlim = c(-3, 3), main = "Dati osservati")
  temp = rmsMK(nsim, mcmc.msmk$theta[[nsim]], mcmc.msmk$sigma[[nsim]], mcmc.msmk$thrs, mcmc.msmk$prob[[nsim]])
  smoothScatter(temp, ylim = c(-3, 3), xlim = c(-3, 3), main = "Dati generati")
  par(mfrow = c(1,1))
}

lpml(y.stand, df = df)
msMK.lpml(y.stand, mcmc.msmk)
