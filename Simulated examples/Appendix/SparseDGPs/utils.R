# Functions Utils

BScore <- function(x,y) mean((x-y)^2)
PEHE <- function(x,y) sqrt(mean((x-y)^2))
MC_se <- function(x, B) qt(0.975, B-1)*sd(x)/sqrt(B)


get_features <- function(N, P) {
  
  # Generate correlated uniforms from a Gaussian Copula
  mysigma = matrix(1, P, P)
  
  for (i in 1:P) {
    for (j in 1:P) {
      mysigma[i, j] = 0.3^abs(i - j) + ifelse(i == j, 0, 0.1)
    }
  }
  
  mycop = MASS::mvrnorm(N, rep(0, P), Sigma = mysigma)
  unif = pnorm(mycop)
  
  # Transform in continuous and binary covariates
  X = matrix(NA, N, P)
  X[, 1:ceiling(P/2)] = qnorm(unif[, 1:ceiling(P/2)])
  X[, (ceiling(P/2)+1):P] = qbinom(unif[, (ceiling(P/2)+1):P], 1, 0.3)
  
  return(X)
  
}


# No Sparse
NOsp_dgp <- function(N, X) {

  # Simulate Pscore and Z
  und_lin = -0.2 + 0.8*X[, 1] - 0.1*X[, 2] + 0.1*X[, 3]*X[, 4] - 0.4*X[, 5] + runif(N)/10
  pscore <- pnorm( und_lin )
  
  Z = rbinom(N, 1, pscore)
  
  ####### RESPONSE SURFACES
  mu = 3 + 1.5*sin(pi*X[, 1]) + 0.5*(X[, 2] - 0.5)^2 + 1.5*(2-abs(X[, 3])) + 1.5*X[, 4]*(X[, 5] + 1)
  ITE = 0.1 +  1*abs(X[, 1] - 1)*(X[, 5] + 2) - 0.4*X[, 3] + 0.6*X[, 2]*X[, 4]
  
  sigma = 1
  Y = mu + ITE*Z + rnorm(N, 0, sigma)
  
  return(list(Y=Y, Z=Z, ITE=ITE))
}


# Sparse PS
sp_PS_dgp <- function(N, X) {
  
  # Simulate Pscore and Z
  und_lin = -0.2 + 0.8*X[, 1] + runif(N)/10
  pscore <- pnorm( und_lin )
  
  Z = rbinom(N, 1, pscore)
  
  ####### RESPONSE SURFACES
  mu = 3 + 1.5*sin(pi*X[, 1]) + 0.5*(X[, 2] - 0.5)^2 + 1.5*(2-abs(X[, 3])) + 1.5*X[, 4]*(X[, 5] + 1)
  ITE = 0.1 +  1*abs(X[, 1] - 1)*(X[, 5] + 2) - 0.4*X[, 3] + 0.6*X[, 2]*X[, 4]
  
  sigma = 1
  Y = mu + ITE*Z + rnorm(N, 0, sigma)
  
  return(list(Y=Y, Z=Z, ITE=ITE))
}


# Sparse Mu
sp_Mu_dgp <- function(N, X) {
  
  # Simulate Pscore and Z
  und_lin = -0.2 + 0.8*X[, 1] - 0.1*X[, 2] + 0.1*X[, 3]*X[, 4] - 0.4*X[, 5] + runif(N)/10
  pscore <- pnorm( und_lin )
  
  Z = rbinom(N, 1, pscore)
  
  ####### RESPONSE SURFACES
  mu = 3 + 1.5*(2-abs(X[, 3]))
  ITE = 0.1 +  1*abs(X[, 1] - 1)*(X[, 5] + 2) - 0.4*X[, 3] + 0.6*X[, 2]*X[, 4]
  
  sigma = 1
  Y = mu + ITE*Z + rnorm(N, 0, sigma)
  
  return(list(Y=Y, Z=Z, ITE=ITE))
}


# Sparse Tau
sp_Tau_dgp <- function(N, X) {
  
  # Simulate Pscore and Z
  und_lin = -0.2 + 0.8*X[, 1] - 0.1*X[, 2] + 0.1*X[, 3]*X[, 4] - 0.4*X[, 5] + runif(N)/10
  pscore <- pnorm( und_lin )
  
  Z = rbinom(N, 1, pscore)
  
  ####### RESPONSE SURFACES
  mu = 3 + 1.5*sin(pi*X[, 1]) + 0.5*(X[, 2] - 0.5)^2 + 1.5*(2-abs(X[, 3])) + 1.5*X[, 4]*(X[, 5] + 1)
  ITE = 0.1 +  1*abs(X[, 1] - 1)
  
  sigma = 1
  Y = mu + ITE*Z + rnorm(N, 0, sigma)
  
  return(list(Y=Y, Z=Z, ITE=ITE))
}
