library("mvnfast")
library("LaplacesDemon")
library("sn")
library("permute")

f1 = function(x, p) dmvn(x, mu = rep(0, p), sigma = diag(1,p))
r1 = function(n, p) rmvn(n, mu = c(0,0), sigma = diag(1,p))


f2 = function(x, p) dmixn(x,
                       rbind(rep(0, p), rep(1/2, p), rep(3/2, p)),
                       list(diag(1, p), diag(2/3, p), diag(5/9, p)),
                       w = (1:p)/(sum(1:p)))
r2 = function(n, p) rmixn(n,
                       rbind(rep(0, p), rep(1/2, p), rep(3/2, p)),
                       list(diag(1, p), diag(2/3, p), diag(5/9, p)),
                       w = (1:p)/(sum(1:p)))


# ---- Mistura 0.5 * N(0, 1/3), 0.5 * N(2, 1/10)
f5_2 = function(x, p=2) dmixn(x,
                       rbind(rep(0, p), rep(1, p)),
                       list(diag(2/3, p), diag(1/10, p)),
                       w = c(5/10, 5/10))
r5_2 = function(n, p=2) rmixn(n,
                       rbind(rep(0, p), rep(1, p)),
                       list(diag(2/3, p), diag(1/10, p)),
                       w = c(5/10, 5/10))

f5_p = function(x, p) dmixn(x,
                       rbind(rep(0, p), rep(1, p)),
                       list(diag(2/3, p), diag(1/10, p)),
                       w = c(5/10, 5/10))
r5_p = function(n, p) rmixn(n,
                       rbind(rep(0, p), rep(1, p)),
                       list(diag(2/3, p), diag(1/10, p)),
                       w = c(5/10, 5/10))

# ---- Mistura 0.5 * SN + 0.5 * N(0, I)
f6_2 = function(x){
    0.5 * dmsn(x,
         xi = rep(0, 2),
         Omega = matrix(0.5, nrow = 2, ncol = 2) + diag(1.5, nrow = 2, ncol = 2),
         alpha = c(-8, 2)) +
    0.5 * dmvn(x, rep(1, 2), diag(1, nrow = 2, ncol = 2))
}

r6_2 = function(n){
  idx = rbinom(n, 1, 0.5)
  n1 = sum(idx == 1)
  return(rbind(
    rmsn(n1,
         xi = rep(0, 2),
         Omega = matrix(0.5, nrow = 2, ncol = 2) + diag(1.5, nrow = 2, ncol = 2),
         alpha = c(-8, 2)
    ),
    rmvn(n - n1, rep(1, 2), diag(1, nrow = 2, ncol = 2)))
  )
}

f6_p = function(x, p = NCOL(x)){
    0.5 * dmsn(x,
         xi = rep(0, p),
         Omega = matrix(0.5, nrow = p, ncol = p) + diag(1.5, nrow = p, ncol = p),
         alpha = c(-8, 2, rep(0, p-2))) +
    0.5 * dmvn(x, rep(1, p), diag(1, nrow = p, ncol = p))
}

r6_p = function(n, p = 2){
  idx = rbinom(n, 1, 0.5)
  n1 = sum(idx == 1)
  return(rbind(
    rmsn(n1,
         xi = rep(0, p),
         Omega = matrix(0.5, nrow = p, ncol = p) + diag(1.5, nrow = p, ncol = p),
         alpha = c(-8, 2, rep(0, p-2))
    ),
    rmvn(n - n1, rep(1, p), diag(1, nrow = p, ncol = p)))
  )
}

r7_2 = function(n){
  idx = rbinom(n, 2, 0.4)
  n2 = sum(idx == 2)
  n1 = sum(idx == 1)
  n0 = sum(idx == 0)
  out = rbind(
    rmsn(n0,
         xi = rep(-2, 2),
         Omega = matrix(1, nrow = 2, ncol = 2) + diag(3, nrow = 2, ncol = 2),
         alpha = c(1, 2)
    ),
    rmsn(n1,
         xi = rep(0, 2),
         Omega = matrix(0.5, nrow = 2, ncol = 2) + diag(1, nrow = 2, ncol = 2),
         alpha = c(-5, 1)
    ),
    rmvn(n2, rep(2, 2), diag(2, nrow = 2, ncol = 2)))
  
  return(out[sample(1:NROW(out), NROW(out)), ])
}

f7_2 = function(x){
  p = dbinom(0:2, 2, 0.4)
  out = p[1] * dmsn(x,
                    xi = rep(-2, 2),
                    Omega = matrix(1, nrow = 2, ncol = 2) + diag(3, nrow = 2, ncol = 2),
                    alpha = c(1, 2)) +
    p[2] * dmsn(x,
                xi = rep(0, 2),
                Omega = matrix(0.5, nrow = 2, ncol = 2) + diag(1, nrow = 2, ncol = 2),
                alpha = c(-5, 1)
                ) +
    p[3] * dmvn(x, rep(2, 2), diag(2, nrow = 2, ncol = 2))
  
  return(out)
}

r7_p = function(n, p=2){
  idx = rbinom(n, 2, 0.4)
  n2 = sum(idx == 2)
  n1 = sum(idx == 1)
  n0 = sum(idx == 0)
  out = rbind(
    rmsn(n0,
         xi = rep(-2, p),
         Omega = matrix(1, nrow = p, ncol = p) + diag(3, nrow = p, ncol = p),
         alpha = c(1, 2, rep(0, p-2))
    ),
    rmsn(n1,
         xi = rep(0, p),
         Omega = matrix(0.5, nrow = p, ncol = p) + diag(1, nrow = p, ncol = p),
         alpha = c(-5, 1, rep(0, p-2))
    ),
    rmvn(n2, rep(2, p), diag(2, nrow = p, ncol = p))
  )
  
  return(out[sample(1:NROW(out), NROW(out)), ])
}

f7_p = function(x, p=NCOL(x)){
  pi = dbinom(0:2, 2, 0.4)
  out = pi[1] * dmsn(x,
                    xi = rep(-2, p),
                    Omega = matrix(1, nrow = p, ncol = p) + diag(3, nrow = p, ncol = p),
                    alpha = c(1, 2, rep(0, p-2))
                    ) +
    pi[2] * dmsn(x,
                xi = rep(0, p),
                Omega = matrix(0.5, nrow = p, ncol = p) + diag(1, nrow = p, ncol = p),
                alpha = c(-5, 1, rep(0, p-2))
                ) +
    pi[3] * dmvn(x, rep(2, p), diag(2, nrow = p, ncol = p))
  
  return(out)
}

r8_2 = function(n){
  idx = rbinom(n, 3, 0.5)
  n3 = sum(idx == 3)
  n2 = sum(idx == 2)
  n1 = sum(idx == 1)
  n0 = sum(idx == 0)
  
  out = rbind(
    rmvn(n0, rep(0, 2), diag(1, nrow = 2, ncol = 2)),
    rmvn(n1, c(3, -1), matrix(0.2,nrow=2,ncol=2) + diag(0.5, nrow = 2, ncol = 2)),
    rmvn(n2, rep(1.5, 2), matrix(0.5, nrow=2, ncol=2) + diag(1.5, nrow = 2, ncol = 2)),
    rmvn(n3, c(0, 2), matrix(0.1, nrow = 2, ncol = 2) + diag(0.2, nrow=2, ncol=2))
    )
  
  return(out[sample(1:NROW(out), NROW(out)), ])
}

f8_2 = function(x){
  p = dbinom(0:3, 3, 0.5)
  out = p[1] * dmvn(x, rep(0, 2), diag(1, nrow = 2, ncol = 2)) +
    p[2] * dmvn(x, c(3, -1), matrix(0.2,nrow=2,ncol=2) + diag(0.5, nrow = 2, ncol = 2)) +
    p[3] * dmvn(x, rep(1.5, 2), matrix(0.5, nrow=2, ncol=2) + diag(1.5, nrow = 2, ncol = 2)) +
    p[4] * dmvn(x, c(0, 2), matrix(0.1, nrow = 2, ncol = 2) + diag(0.2, nrow=2, ncol=2))
  return(out)
}

r8_p = function(n, p=2){
  idx = rbinom(n, 3, 0.5)
  n3 = sum(idx == 3)
  n2 = sum(idx == 2)
  n1 = sum(idx == 1)
  n0 = sum(idx == 0)
  
  out = rbind(
    rmvn(n0, rep(0, p), diag(1, nrow = p, ncol = p)),
    rmvn(n1, c(3, -1, rep(0, p-2)), matrix(0.2,nrow=p,ncol=p) + diag(0.5, nrow = p, ncol = p)),
    rmvn(n2, rep(1.5, p), matrix(0.5, nrow=p, ncol=p) + diag(1.5, nrow = p, ncol = p)),
    rmvn(n3, c(0, 2, rep(0, p-2)), matrix(0.1, nrow = p, ncol = p) + diag(0.2, nrow=p, ncol=p))
    )
  
  return(out[sample(1:NROW(out), NROW(out)), ])
}

f8_p = function(x, p=NCOL(x)){
  pi = dbinom(0:3, 3, 0.5)
  out = pi[1] * dmvn(x, rep(0, p), diag(1, nrow = p, ncol = p)) +
    pi[2] * dmvn(x, c(3, -1, rep(0, p-2)), matrix(0.2,nrow=p,ncol=p) + diag(0.5, nrow = p, ncol = p)) +
    pi[3] * dmvn(x, rep(1.5, p), matrix(0.5, nrow=p, ncol=p) + diag(1.5, nrow = p, ncol = p)) +
    pi[4] * dmvn(x, c(0, 2, rep(0, p-2)), matrix(0.1, nrow = p, ncol = p) + diag(0.2, nrow=p, ncol=p))
  return(out)
}

f9_2 = function(x){
  p = c(0.2, 0.3, 0.3, 0.2)
  out = p[1] * dmvn(x, c(-2, -1), diag(0.5, nrow = 2, ncol = 2)) +
    p[2] * dmvn(x, c(0, 1), matrix(0.2, nrow=2,ncol=2) + diag(0.5, nrow = 2, ncol = 2)) +
    p[3] * dmvn(x, c(2, 2), matrix(-0.2, nrow=2, ncol=2) + diag(1, nrow = 2, ncol = 2)) +
    p[4] * dmvn(x, c(1, -1), matrix(0.1, nrow = 2, ncol = 2) + diag(0.2, nrow=2, ncol=2))
  return(out)
}

r9_2 = function(n){
  idx = sample(0:3, size = n, replace = TRUE, prob = c(0.2, 0.3, 0.3, 0.2))
  n3 = sum(idx == 3)
  n2 = sum(idx == 2)
  n1 = sum(idx == 1)
  n0 = sum(idx == 0)
  
  out = rbind(
    rmvn(n0, c(-2, -1), diag(0.5, nrow = 2, ncol = 2)),
    rmvn(n1, c(0, 1), matrix(0.2, nrow=2,ncol=2) + diag(0.5, nrow = 2, ncol = 2)),
    rmvn(n2, c(2, 2), matrix(-0.2, nrow=2, ncol=2) + diag(1, nrow = 2, ncol = 2)),
    rmvn(n3, c(1, -1), matrix(0.1, nrow = 2, ncol = 2) + diag(0.2, nrow=2, ncol=2))
    )
  
  return(out[sample(1:NROW(out), NROW(out)), ])
}

f9_p = function(x, p = NCOL(x)){
  pi = c(0.2, 0.3, 0.3, 0.2)
  out = pi[1] * dmvn(x, c(-2, -1, rep(0, p - 2)), diag(0.5, nrow = p, ncol = p)) +
    pi[2] * dmvn(x, c(0, 1, rep(0, p-2)), matrix(0.2, nrow=p,ncol=p) + diag(0.5, nrow = p, ncol = p)) +
    pi[3] * dmvn(x, c(2, 2, rep(0, p-2)), matrix(0.2, nrow=p, ncol=p) + diag(1, nrow = p, ncol = p)) +
    pi[4] * dmvn(x, c(1, -1, rep(0, p-2)), matrix(0.1, nrow = p, ncol = p) + diag(0.2, nrow=p, ncol=p))
  return(out)
}

r9_p = function(n, p = 2){
  idx = sample(0:3, size = n, replace = TRUE, prob = c(0.2, 0.3, 0.3, 0.2))
  n3 = sum(idx == 3)
  n2 = sum(idx == 2)
  n1 = sum(idx == 1)
  n0 = sum(idx == 0)
  
  out = rbind(
    rmvn(n0, c(-2, -1, rep(0, p-2)), diag(0.5, nrow = p, ncol = p)),
    rmvn(n1, c(0, 1, rep(0, p-2)), matrix(0.2, nrow=p,ncol=p) + diag(0.5, nrow = p, ncol = p)),
    rmvn(n2, c(2, 2, rep(0, p-2)), matrix(0.2, nrow=p, ncol=p) + diag(1, nrow = p, ncol = p)),
    rmvn(n3, c(1, -1, rep(0, p-2)), matrix(0.1, nrow = p, ncol = p) + diag(0.2, nrow=p, ncol=p))
    )
  
  return(out[sample(1:NROW(out), NROW(out)), ])
}

f10_2 = function(x){
  p = c(4, 1, 1, 1, 1, 1, 2, 1)
  p = p/sum(p)
  out = p[1] * dmvn(x, c(0, 0), diag(1, nrow = 2, ncol = 2)) +
    p[2] * dmvn(x, c(1, 1), matrix(0.015, nrow=2,ncol=2) + diag(0.05, nrow = 2, ncol = 2)) +
    p[3] * dmvn(x, c(1, -1), matrix(-0.01, nrow=2, ncol=2) + diag(0.05, nrow = 2, ncol = 2)) +
    p[4] * dmvn(x, c(0, -2), matrix(0.05, nrow = 2, ncol = 2) + diag(0.1, nrow=2, ncol=2)) +
    p[5] * dmvn(x, c(-1, 2), matrix(0.1, nrow = 2, ncol = 2) + diag(0.25, nrow=2, ncol=2)) +
    p[6] * dmvn(x, c(-1, -2), matrix(-0.05, nrow = 2, ncol = 2) + diag(0.25, nrow=2, ncol=2)) +
    p[7] * dmvn(x, c(-1.5, 0), matrix(0, nrow = 2, ncol = 2) + diag(0.25, nrow=2, ncol=2)) +
    p[8] * dmsn(x, xi = c(0,0), Omega = diag(0.25, nrow=2, ncol=2), alpha = c(5, 4))
    
  return(out)
}

r10_2 = function(n){
  p = c(4, 1, 1, 1, 1, 1, 2, 1)
  p = p/sum(p)
  idx = sample(0:(length(p)-1), size = n, replace = TRUE, prob = p)
  n7 = sum(idx == 7)
  n6 = sum(idx == 6)
  n5 = sum(idx == 5)
  n4 = sum(idx == 4)
  n3 = sum(idx == 3)
  n2 = sum(idx == 2)
  n1 = sum(idx == 1)
  n0 = sum(idx == 0)
  
  out = rbind(
    rmvn(n0, c(0, 0), diag(0.75, nrow = 2, ncol = 2)) ,
    rmvn(n1, c(1, 1), matrix(0.015, nrow=2,ncol=2) + diag(0.05, nrow = 2, ncol = 2)), 
    rmvn(n2, c(1, -1), matrix(-0.01, nrow=2, ncol=2) + diag(0.05, nrow = 2, ncol = 2)) ,
    rmvn(n3, c(0, -2), matrix(0.05, nrow = 2, ncol = 2) + diag(0.1, nrow=2, ncol=2)) ,
    rmvn(n4, c(-1, 2), matrix(0.1, nrow = 2, ncol = 2) + diag(0.25, nrow=2, ncol=2)) ,
    rmvn(n5, c(-1, -2), matrix(-0.05, nrow = 2, ncol = 2) + diag(0.25, nrow=2, ncol=2)), 
    rmvn(n6, c(-1.5, 0), matrix(0, nrow = 2, ncol = 2) + diag(0.25, nrow=2, ncol=2)),
    rmsn(n7, xi = c(0,0), Omega = diag(0.25, nrow=2, ncol=2), alpha = c(5, 4))
    )
  
  return(out[sample(1:NROW(out), NROW(out)), ])
}

f10_p = function(x, p = NCOL(x)){
  pi = c(4, 1, 1, 1, 1, 1, 2, 1)
  pi = pi/sum(pi)
  out = pi[1] * dmvn(x, c(0, 0, rep(0, p-2)), diag(1, nrow = p, ncol = p)) +
    pi[2] * dmvn(x, c(1, 1, rep(0, p-2)), matrix(0.015, nrow=p,ncol=p) + diag(0.15, nrow = p, ncol = p)) +
    pi[3] * dmvn(x, c(1, -1, rep(0, p-2)), matrix(-0.01, nrow=p, ncol=p) + diag(0.15, nrow = p, ncol = p)) +
    pi[4] * dmvn(x, c(0, -2, rep(0, p-2)), matrix(0.05, nrow = p, ncol = p) + diag(0.1, nrow=p, ncol=p)) +
    pi[5] * dmvn(x, c(-1, 2, rep(0, p-2)), matrix(0.1, nrow = p, ncol = p) + diag(0.25, nrow=p, ncol=p)) +
    pi[6] * dmvn(x, c(-1, -2, rep(0, p-2)), matrix(0.05, nrow = p, ncol = p) + diag(0.25, nrow=p, ncol=p)) +
    pi[7] * dmvn(x, c(-1.5, 0, rep(0, p-2)), matrix(0, nrow = p, ncol = p) + diag(0.25, nrow=p, ncol=p)) +
    pi[8] * dmsn(x, xi = c(0,0, rep(0, p-2)), Omega = diag(0.25, nrow=p, ncol=p), alpha = c(5, 4, rep(0, p-2)))
    
  return(out)
}

r10_p = function(n, p){
  pi = c(4, 1, 1, 1, 1, 1, 2, 1)
  pi = pi/sum(pi)
  idx = sample(0:(length(pi)-1), size = n, replace = TRUE, prob = pi)
  n7 = sum(idx == 7)
  n6 = sum(idx == 6)
  n5 = sum(idx == 5)
  n4 = sum(idx == 4)
  n3 = sum(idx == 3)
  n2 = sum(idx == 2)
  n1 = sum(idx == 1)
  n0 = sum(idx == 0)
  
  out = rbind(
    rmvn(n0, c(0, 0, rep(0, p-2)), diag(0.75, nrow = p, ncol = p)) ,
    rmvn(n1, c(1, 1, rep(0, p-2)), matrix(0.015, nrow=p,ncol=p) + diag(0.15, nrow = p, ncol = p)), 
    rmvn(n2, c(1, -1, rep(0, p-2)), matrix(-0.01, nrow=p, ncol=p) + diag(0.15, nrow = p, ncol = p)) ,
    rmvn(n3, c(0, -2, rep(0, p-2)), matrix(0.05, nrow = p, ncol = p) + diag(0.1, nrow=p, ncol=p)) ,
    rmvn(n4, c(-1, 2, rep(0, p-2)), matrix(0.1, nrow = p, ncol = p) + diag(0.25, nrow=p, ncol=p)) ,
    rmvn(n5, c(-1, -2, rep(0, p-2)), matrix(0.05, nrow = p, ncol = p) + diag(0.25, nrow=p, ncol=p)), 
    rmvn(n6, c(-1.5, 0, rep(0, p-2)), matrix(0, nrow = p, ncol = p) + diag(0.25, nrow=p, ncol=p)),
    rmsn(n7, xi = c(0,0, rep(0, p-2)), Omega = diag(0.25, nrow=p, ncol=p), alpha = c(5, 4, rep(0, p-2)))
    )
  
  return(out[sample(1:NROW(out), NROW(out)), ])
}

f11_2 = function(x){
  p = c(20, rep(1, 7))
  p = p/sum(p)
  out = p[1] * mvnfast::dmvn(x, c(0, 0), diag(2.5, nrow = 2, ncol = 2) + matrix(0.2, nrow=2, ncol=2)) +
    p[2] * dmvn(x, c(1, 1), matrix(0.015, nrow=2,ncol=2) + diag(0.05, nrow = 2, ncol = 2)) +
    p[3] * dmvn(x, c(1, -1),  diag(0.03, nrow = 2, ncol = 2)) +
    p[4] * dmvn(x, c(-1, 0), matrix(-0.015, nrow=2,ncol=2) + diag(0.05, nrow = 2, ncol = 2)) +
    p[5] * dmvn(x, c(2, 2), matrix(0.035, nrow=2,ncol=2) + diag(0.1, nrow = 2, ncol = 2)) +
    p[6] * dmvn(x, c(2, -2), matrix(-0.035, nrow=2,ncol=2) + diag(0.1, nrow = 2, ncol = 2)) +
    p[7] * dmvn(x, c(-2, 2), matrix(-0.035, nrow=2,ncol=2) + diag(0.1, nrow = 2, ncol = 2)) +
    p[8] * dmvn(x, c(-2, -2), matrix(0.035, nrow=2,ncol=2) + diag(0.1, nrow = 2, ncol = 2)) 
    
  return(out)
}

r11_2 = function(n){
  p = c(20, rep(1, 7))
  p = p/sum(p)
  idx = sample(0:(length(p)-1), size = n, replace = TRUE, prob = p)
  n7 = sum(idx == 7)
  n6 = sum(idx == 6)
  n5 = sum(idx == 5)
  n4 = sum(idx == 4)
  n3 = sum(idx == 3)
  n2 = sum(idx == 2)
  n1 = sum(idx == 1)
  n0 = sum(idx == 0)
  
  out = rbind(mvnfast::rmvn(n0, c(0, 0), diag(2.5, nrow = 2, ncol = 2) + matrix(0.2, nrow=2, ncol=2)),
    rmvn(n1, c(1, 1), matrix(0.015, nrow=2,ncol=2) + diag(0.05, nrow = 2, ncol = 2)),
    rmvn(n2, c(1, -1), diag(0.03, nrow = 2, ncol = 2)),
    rmvn(n3, c(-1, 0), matrix(-0.015, nrow=2,ncol=2) + diag(0.05, nrow = 2, ncol = 2)),
    rmvn(n4, c(2, 2), matrix(0.035, nrow=2,ncol=2) + diag(0.1, nrow = 2, ncol = 2)),
    rmvn(n5, c(2, -2), matrix(-0.035, nrow=2,ncol=2) + diag(0.1, nrow = 2, ncol = 2)),
    rmvn(n6, c(-2, 2), matrix(-0.035, nrow=2,ncol=2) + diag(0.1, nrow = 2, ncol = 2)),
    rmvn(n7, c(-2, -2), matrix(0.035, nrow=2,ncol=2) + diag(0.1, nrow = 2, ncol = 2)) 
    )
  
  return(out[sample(1:NROW(out), NROW(out)), ])
}

f11_p = function(x, p = NCOL(x)){
  pi = c(20, rep(1, 7))
  pi = pi/sum(pi)
  out = pi[1] * mvnfast::dmvn(x, c(0, 0, rep(0, p-2)), diag(2.5, nrow = p, ncol = p) + matrix(0.2, nrow=p, ncol=p)) +
    pi[2] * dmvn(x, c(1, 1, rep(0, p-2)), matrix(0.015, nrow=p,ncol=p) + diag(0.05, nrow = p, ncol = p)) +
    pi[3] * dmvn(x, c(1, -1, rep(0, p-2)),  diag(0.03, nrow = p, ncol = p)) +
    pi[4] * dmvn(x, c(-1, 0, rep(0, p-2)), matrix(0.015, nrow=p,ncol=p) + diag(0.05, nrow = p, ncol = p)) +
    pi[5] * dmvn(x, c(2, 2, rep(0, p-2)), matrix(0.035, nrow=p,ncol=p) + diag(0.1, nrow = p, ncol = p)) +
    pi[6] * dmvn(x, c(2, -2, rep(0, p-2)), matrix(0.035, nrow=p,ncol=p) + diag(0.1, nrow = p, ncol = p)) +
    pi[7] * dmvn(x, c(-2, 2, rep(0, p-2)), matrix(0.035, nrow=p,ncol=p) + diag(0.1, nrow = p, ncol = p)) +
    pi[8] * dmvn(x, c(-2, -2, rep(0, p-2)), matrix(0.035, nrow=p,ncol=p) + diag(0.1, nrow = p, ncol = p)) 
    
  return(out)
}

r11_p = function(n, p){
  pi = c(20, rep(1, 7))
  pi = pi/sum(pi)
  idx = sample(0:(length(pi)-1), size = n, replace = TRUE, prob = pi)
  n7 = sum(idx == 7)
  n6 = sum(idx == 6)
  n5 = sum(idx == 5)
  n4 = sum(idx == 4)
  n3 = sum(idx == 3)
  n2 = sum(idx == 2)
  n1 = sum(idx == 1)
  n0 = sum(idx == 0)
  
  out = rbind(mvnfast::rmvn(n0, c(0, 0, rep(0, p-2)), diag(2.5, nrow = p, ncol = p) + matrix(0.2, nrow=p, ncol=p)),
    rmvn(n1, c(1, 1, rep(0, p-2)), matrix(0.015, nrow=p,ncol=p) + diag(0.05, nrow = p, ncol = p)),
    rmvn(n2, c(1, -1, rep(0, p-2)), diag(0.03, nrow = p, ncol = p)),
    rmvn(n3, c(-1, 0, rep(0, p-2)), matrix(0.015, nrow=p,ncol=p) + diag(0.05, nrow = p, ncol = p)),
    rmvn(n4, c(2, 2, rep(0, p-2)), matrix(0.035, nrow=p,ncol=p) + diag(0.1, nrow = p, ncol = p)),
    rmvn(n5, c(2, -2, rep(0, p-2)), matrix(0.035, nrow=p,ncol=p) + diag(0.1, nrow = p, ncol = p)),
    rmvn(n6, c(-2, 2, rep(0, p-2)), matrix(0.035, nrow=p,ncol=p) + diag(0.1, nrow = p, ncol = p)),
    rmvn(n7, c(-2, -2, rep(0, p-2)), matrix(0.035, nrow=p,ncol=p) + diag(0.1, nrow = p, ncol = p)) 
    )
  
  return(out[sample(1:NROW(out), NROW(out)), ])
}

f12_2 = function(x){
  p = c(6, rep(1, 5))
  p = p/sum(p)
  out = p[1] * mvnfast::dmvn(x, c(0, 0), diag(2.5, nrow = 2, ncol = 2) + matrix(0.5, nrow=2, ncol=2)) +
    p[2] * dmvn(x, c(-2, 1), matrix(0.0015, nrow=2,ncol=2) + diag(0.005, nrow = 2, ncol = 2)) +
    p[3] * dmvn(x, c(1, -1),  diag(0.008, nrow = 2, ncol = 2)) +
    p[4] * dmvn(x, c(2, 1.5), matrix(0.015, nrow=2,ncol=2) + diag(0.05, nrow = 2, ncol = 2)) +
    p[5] * dmvn(x, c(2, 2), matrix(0.035, nrow=2,ncol=2) + diag(0.1, nrow = 2, ncol = 2)) +
    p[6] * dmvn(x, c(2, -2), matrix(0.035, nrow=2,ncol=2) + diag(0.1, nrow = 2, ncol = 2))
  return(out)
}

r12_2 = function(n){
  p = c(6, rep(1, 5))
  p = p/sum(p)
  idx = sample(0:(length(p)-1), size = n, replace = TRUE, prob = p)
  n5 = sum(idx == 5)
  n4 = sum(idx == 4)
  n3 = sum(idx == 3)
  n2 = sum(idx == 2)
  n1 = sum(idx == 1)
  n0 = sum(idx == 0)
  
  out = rbind(mvnfast::rmvn(n0, c(0, 0), diag(2.5, nrow = 2, ncol = 2) + matrix(0.5, nrow=2, ncol=2)),
    rmvn(n1, c(-2, 1), matrix(0.0015, nrow=2,ncol=2) + diag(0.005, nrow = 2, ncol = 2)),
    rmvn(n2, c(1, -1), diag(0.008, nrow = 2, ncol = 2)),
    rmvn(n3, c(2, 1.5), matrix(0.015, nrow=2,ncol=2) + diag(0.05, nrow = 2, ncol = 2)),
    rmvn(n4, c(2, 2), matrix(0.035, nrow=2,ncol=2) + diag(0.1, nrow = 2, ncol = 2)),
    rmvn(n5, c(2, -2), matrix(0.035, nrow=2,ncol=2) + diag(0.1, nrow = 2, ncol = 2))
    )
  
  return(out[sample(1:NROW(out), NROW(out)), ])
}

f12_p = function(x, p){
  pi = c(6, rep(1, 5))
  pi = pi/sum(pi)
  out = pi[1] * mvnfast::dmvn(x, c(0, 0, rep(0, p-2)), diag(2.5, nrow = p, ncol = p) + matrix(0.5, nrow=p, ncol=p)) +
    pi[2] * dmvn(x, c(-2, 1, rep(0, p-2)), matrix(0.0015, nrow=p,ncol=p) + diag(0.005, nrow = p, ncol = p)) +
    pi[3] * dmvn(x, c(1, -1, rep(0, p-2)),  diag(0.008, nrow = p, ncol = p)) +
    pi[4] * dmvn(x, c(2, 1.5, rep(0, p-2)), matrix(0.015, nrow=p,ncol=p) + diag(0.05, nrow = p, ncol = p)) +
    pi[5] * dmvn(x, c(2, 2, rep(0, p-2)), matrix(0.035, nrow=p,ncol=p) + diag(0.1, nrow = p, ncol = p)) +
    pi[6] * dmvn(x, c(2, -2, rep(0, p-2)), matrix(0.035, nrow=p,ncol=p) + diag(0.1, nrow = p, ncol = p))
  return(out)
}

r12_p = function(n, p){
  pi = c(6, rep(1, 5))
  pi = pi/sum(pi)
  idx = sample(0:(length(pi)-1), size = n, replace = TRUE, prob = pi)
  n5 = sum(idx == 5)
  n4 = sum(idx == 4)
  n3 = sum(idx == 3)
  n2 = sum(idx == 2)
  n1 = sum(idx == 1)
  n0 = sum(idx == 0)
  
  out = rbind(mvnfast::rmvn(n0, c(0, 0, rep(0, p-2)), diag(2.5, nrow = p, ncol = p) + matrix(0.5, nrow=p, ncol=p)),
    rmvn(n1, c(-2, 1, rep(0, p-2)), matrix(0.0015, nrow=p,ncol=p) + diag(0.005, nrow = p, ncol = p)),
    rmvn(n2, c(1, -1, rep(0, p-2)), diag(0.008, nrow = p, ncol = p)),
    rmvn(n3, c(2, 1.5, rep(0, p-2)), matrix(0.015, nrow=p,ncol=p) + diag(0.05, nrow = p, ncol = p)),
    rmvn(n4, c(2, 2, rep(0, p-2)), matrix(0.035, nrow=p,ncol=p) + diag(0.1, nrow = p, ncol = p)),
    rmvn(n5, c(2, -2, rep(0, p-2)), matrix(0.035, nrow=p,ncol=p) + diag(0.1, nrow = p, ncol = p))
    )
  
  return(out[sample(1:NROW(out), NROW(out)), ])
}
#---- Density contour for a function f(x)
densContour = function(f, y = NULL, npoints = 100,  xlim = c(-3, 3), ylim = c(-3, 3), ...){
  if(!is.null(y)){
    xx = seq(min(y[ , 1]), max(y[ , 1]), length.out = npoints)
    yy = seq(min(y[ , 2]), max(y[ , 2]), length.out = npoints)
  } else{
    xx = seq(xlim[1], xlim[2], length.out = npoints)
    yy = seq(ylim[1], ylim[2], length.out = npoints)
  }
    grid = as.matrix(expand.grid(xx,yy))
    zz = matrix(f(grid), nrow = npoints)
    contour(xx, yy, zz, ...)
}


#---- Calculate alpha from delta for expected depth
deltaToAlpha = function(delta){
  if(delta == 0) return(c(1, 3, 5))
  if(delta == 0.25) return(c(0.25, 1.25, 2.25))
  if(delta == 0.5) return(c(-0.45, -0.35, -0.25))
}

