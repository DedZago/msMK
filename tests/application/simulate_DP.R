library(DPpackage)
setwd("/home/dede/Documents/git/MSc-thesis/code/tests/application")
nsim = 1000
burnin = 1000
dati = read.csv("new_galaxy_data.csv")

y = dati[ , -5 ]
y.stand = y
y.stand = as.matrix(y.stand)
  
# ---- Create splits
p = NCOL(y.stand)
b = 1
k0 = apply(y, 2, var)

# Initial state
state <- NULL

# MCMC parameters
mu0 = as.vector(apply(y, 2, mean))
sig0 = diag(apply(y, 2, var))

nburn <- burnin
nsave <- nsim
nskip <- 0
ndisplay <- 100
mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)

prior <- list(alpha=2.5, m1=mu0, psiinv1=sig0,
               nu1=p, k0 = 1)
#DP plus hyperprior
# prior <- list(a0=2,b0=2, m2=mu0, psiinv2=sig0,
#               nu1=p+1, nu2=p+1, 
#               s2 = diag(2, nrow=p),
#               tau1=1,tau2=100)

fit1.1 <- NULL
attempt <- 1
while( is.null(fit1.1) && attempt <= 50 ) {
  cat("Attempt", attempt, "\n")
  attempt <- attempt + 1
  try(
      fit1.1 <- DPdensity(y.stand, prior=prior,mcmc=mcmc, state=state, status=TRUE)
  )
} 

cat("Saving model image\n\n")
# save.image("model_estimated_DP.Rdata")



state = fit1.1$state
fit1.1 <- DPdensity(y.stand, prior=prior,mcmc=mcmc, state=state, status=FALSE)
load("model_estimated_DP.Rdata")
