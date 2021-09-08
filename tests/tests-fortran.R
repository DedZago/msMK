library("devtools")
setwd("~/Documents/git/MSc-thesis/code/msMK")
load_all()
set.seed(194289)

n = 4000
d = 2
mean = c(0,0)
sigma = matrix(c(1, 0, 0, 1), byrow = TRUE, ncol = d)
lower = c(-Inf, -Inf)
upper = c(Inf, Inf)
x0 = c(0,0)
burnin = 1
thinning = 1
# y = rtmvnorm_rcpp(n, mean, sigma, lower, upper, burnin, thinning)
y = rtmvnorm_arma(n, mean, sigma, lower, upper, burnin, thinning)
print(head(y))

y1 = y[ , 1]
y2 = y[ , 2]
##  Create cuts:
x_c <- cut(y[ , 1], 50)
y_c <- cut(y[ , 2], 50)

##  Calculate joint counts at cut levels:
z <- table(x_c, y_c)

library(plot3D)
hist3D(z=z, border="black")
