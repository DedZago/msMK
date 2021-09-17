# Test normal distribution
library(devtools)
library(Rcpp)
library(RcppArmadillo)
# devtools::install_github("lbelzile/TruncatedNormal")
library(TruncatedNormal)
library(tmvtnorm)
library(microbenchmark)
library(kableExtra)
library(ggplot2)
setwd("~/Documents/git/MSc-thesis/code/msMK")
source("tests/tests.R")
setwd("~/Documents/git/MSc-thesis/code/tests/time")

# Wiener process plot
set.seed(123)
T = 1
N = 501
dt = T/N

dW = sqrt(dt)*rnorm(N)   # increments
dX = sqrt(dt)*rnorm(N)
sp = cbind(cumsum(dW), cumsum(dX))           # cumulative sum
matplot(x = seq(0,T, length=NROW(sp)), sp, type="l", lty = 1, col = c("blue", "darkorange"),
        ylab = "X", xlab = "t", lwd=2)




alb = c(3, -2)
ub = c(Inf, Inf)
mu = c(0, 0)
Sigma = matrix(0.5, nrow=2, ncol=2) + diag(3, nrow = 2, ncol=2)

# Botev vs Gibbs sampler from Horrace + Wilhelm
microbenchmark(
  "TruncatedNormal" = TruncatedNormal::rtmvnorm(1000, mu, Sigma, lb, ub),
  "tmvtnorm" = tmvtnorm::rtmvnorm(1000, mu, Sigma, lb, ub, algorithm = "gibbs"),
  times = 1e03
  ) %>%
  summary() %>%
  kbl(format = "latex", booktabs = TRUE, digits = 2) %>%
  kable_styling(latex_options = "HOLD_position")

# Rtruncnorm d singles vs Botev multivariate distribution
d = 2
mu = rep(0, d)
mvar = rep(2, d)
lb = rep(-Inf, d)
ub = rep(Inf, d)
Sigma = diag(mvar, ncol = d, nrow = d)


rBotev = function(s){
   for(i in 1:(2^(s+1)-1)){
     TruncatedNormal::rtmvnorm(1, mu, Sigma, lb, ub)
   }
}

rSingles = function(s){
  for(i in 1:(2^(s+1)-1)){
    for(j in 1:d){
      truncnorm::rtruncnorm(1, mu[j], mvar[j], lb[j], ub[j])
    }
  }
}

smax = 10
times = 30
out = data.frame("s" = NULL, "Singles" = NULL, "Botev" = NULL)
for(s in 1:smax){
  temp =  microbenchmark(
    "Indep" = rSingles(s),
    "Full" = rBotev(s),
    times = times
  )
  medt = quantile(temp$time[temp$expr == "Indep"], probs = c(0.25, 0.5, 0.75)) / 1000
  medd = quantile(temp$time[temp$expr == "Full"], probs = c(0.25, 0.5, 0.75)) / 1000
  out = rbind(out, c(s, medt, medd))
}

colnames(out) = c("s", "Indep.low", "Indep", "Indep.hi", "Full.low", "Full", "Full.hi")
out

ggplot(out, aes(x = s)) +
  geom_point(aes(y = Full, colour = "darkorange")) +
  geom_errorbar(aes(ymin = Full.low, ymax = Full.hi, colour = "darkorange"),
                width=.2, position=position_dodge(0.05)) +
  geom_point(aes(y = Indep, , colour = "blue")) +
  geom_errorbar(aes(ymin = Indep.low, ymax = Indep.hi, colour = "blue"),
                width=.2, position=position_dodge(0.05)) +
  xlab("Maximum depth") +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  ylab("time (ms)") +
  scale_colour_identity(guide = "legend",
                        name = "Base measure",
                        labels = c("Independent", "Full covariance"))+
  theme_bw()
  
ggsave("indepFullComparison.png", device = "png", width=6, height=4, dpi = 600)

# Botev diag vs Botev non diag

rdiagsingle = function(d){
  out = rep(NA, d)
  for(i in 1:d){
    out[i] = TruncatedNormal::rtmvnorm(1, 0, 1, 0, Inf)
  }
  out
}
rdiag = function(d) TruncatedNormal::rtmvnorm(1, rep(0, d), diag(1, nrow=d, ncol=d), rep(0, d), rep(Inf, d))
rfull = function(d) TruncatedNormal::rtmvnorm(1, rep(0, d), matrix(0.5, nrow=d, ncol=d) + diag(1, nrow=d, ncol=d),
                                              rep(0, d), rep(Inf, d))

d = 2
microbenchmark::microbenchmark(
               "DiagSingle" = rdiagsingle(d),
               "Diag" = rdiag(d),
               "Full" = rfull(d),
               times = 100
)

setwd("~/Documents/git/MSc-thesis/code/tests/time")
sourceCpp("comparison.cpp")

smax = 11
times = 30
out = data.frame("s" = NULL, "Tony" = NULL, "Dede" = NULL)
for(s in 1:smax){
  temp =  microbenchmark(
    "Tony" = bintree_tony(s),
    "Dede" = bintree_dede(s),
    times = times
  )
  medt = quantile(temp$time[temp$expr == "Tony"], probs = c(0.25, 0.5, 0.75)) / 1000
  medd = quantile(temp$time[temp$expr == "Dede"], probs = c(0.25, 0.5, 0.75)) / 1000
  out = rbind(out, c(s, medt, medd))
}

save(out, file = "out.Rdata")
colnames(out) = c("s", "Tree.low", "Tree", "Tree.hi", "Array.low", "Array", "Array.hi")
out

library(ggplot2)
load("out.Rdata")
colnames(out) = c("s", "Tree.low", "Tree", "Tree.hi", "Array.low", "Array", "Array.hi")
out

ggplot(out, aes(x = s)) +
  geom_point(aes(y = Tree, colour = "darkorange")) +
  geom_errorbar(aes(ymin = Tree.low, ymax = Tree.hi, colour = "darkorange"),
                width=.2, position=position_dodge(0.05)) +
  geom_point(aes(y = Array, , colour = "blue")) +
  geom_errorbar(aes(ymin = Array.low, ymax = Array.hi, colour = "blue"),
                width=.2, position=position_dodge(0.05)) +
  xlab("Maximum depth") +
  ylab("time (ms)") +
  scale_colour_identity(guide = "legend",
                        name = "Method",
                        labels = c("array", "msBP")) +
  theme_bw()
  
ggsave("tonyDedeComparison.png", device = "png", width=6, height=4, dpi = 600)
  

# msMK indep vs non indep

library("devtools")
library(RcppArmadillo)
library(Rcpp)
library("reticulate")
setwd("~/Documents/git/MSc-thesis/code/msMK")
source("tests/tests.R")
devtools::load_all()

reticulate::py_install(c("numpy", "numpy-hilbert-curve"), pip = TRUE)
reticulate::source_python("hilbertSplit.py")

setwd("~/Documents/git/MSc-thesis/code/tests/time")
delta=0.25
a = deltaToAlpha(delta)[2]
y = r9_2(100)
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
nsim = 5
burnin = 5

microbenchmark::microbenchmark(
  "indep" = indep(),
  "full" = full(),
  times = 10
)

single = function(n, p) for(j in 1:p) truncnorm::rtruncnorm(n, 1, 2, 0, 1)
multi = function(n, p) TruncatedNormal::rtmvnorm(n, rep(0, p), diag(1, nrow=p, ncol=p), rep(1, p), rep(2, p))

microbenchmark::microbenchmark(
  "single" = single(100, 20),
  "multi" = multi(100, 20)
)

dsingle = function(n, p){
  for(j in 1:p) truncnorm::dtruncnorm(rnorm(n), 1, 2, 0, 1)
}

dsingleBotev = function(n, p){
  for(j in 1:p) TruncatedNormal::dtmvnorm(matrix(rnorm(n), ncol=1), 0, 1, 1, 2)
}

dmultiBotev = function(n, p) TruncatedNormal::dtmvnorm(matrix(rnorm(p*n), nrow=n), rep(0, p), diag(1, nrow=p, ncol=p), rep(1, p), rep(2, p))

# Se n > 2900 usare multibotev, altrimenti usare singlebotev
# Forse
p = 2
n = 3000
microbenchmark::microbenchmark(
  "single" = dsingle(n,p),
  "singleBotev" = dsingleBotev(n,p),
  "multiBotev" = dmultiBotev(n,p),
  times = 10
)
sourceCpp("comparison.cpp")

Rfor = function(n){
  a = 0
  for(i in 1:n){
    a = a + 1
  }
}

microbenchmark::microbenchmark(
  "R" = Rfor(10000),
  "C++" = lol(10000),
  times = 10
)



library(ggplot2)
library(ggpubr)
# Dirichlet process examples
rDP = function(n, a, G_0){
  b = rbeta(n, 1, a)
  p = numeric(n)
  p[1] <- b[1]
  p[2:n] <- sapply(2:n, function(i) b[i] * prod(1 - b[1:(i-1)]))
  y = G_0(n)
  return(sample(y, prob = p, replace = TRUE))
}

n = 10000
G_0 = function(n) rnorm(n, 0, 1)

a = c(2, 10, 50)
plots = vector(mode="list", length=length(a))
set.seed(1)
for(i in 1:length(a)){
  df = data.frame("theta" = rDP(n, a[i], G_0))
  df %>%
    ggplot(aes(x = theta)) +
    geom_histogram(bins = 100, col = "white", freq = FALSE) +
    theme_bw() -> plots[[i]] 
}
ggarrange(plotlist=plots, ncol = 1)
ggsave("dpsamples.png", device = "png", width=6.5, height=4.2, dpi = 600)

theta = rDP(n, a, G_0)
hist(theta, nclass = 100, prob = TRUE)

