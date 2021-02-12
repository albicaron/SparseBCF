######################################
# 1st Simulated Example in the Paper #
######################################

rm(list = ls())

# Libraries
library(tidyverse)
library(BART)
library(SparseBCF)

# Functions
bias <- function(x,y) mean(x-y)
PEHE <- function(x,y) sqrt(mean((x-y)^2))
coverage_95 <- function(ITE_EST, ITE) {
  quan = apply(ITE_EST, 2, function(x) quantile(x, c(0.025, 0.975)))
  cove = sum((quan[1,] < ITE) & (ITE < quan[2,]))/length(ITE)
  
  return(cove)
}
MC_se <- function(x, B) qt(0.975, B-1)*sd(x)/sqrt(B)

# Options
P = 10
N = 1000
B = 500


# Store Metrics
Bias_BCF = c(NA); Bias_5BCF = c(NA); Bias_SBCF = c(NA)
PEHE_BCF = c(NA); PEHE_5BCF = c(NA); PEHE_SBCF = c(NA)
Cover_BCF = c(NA); Cover_5BCF = c(NA); Cover_SBCF = c(NA)

Split_Mu = matrix(NA, B, P+1); Split_Tau = matrix(NA, B, P)


for (b in 1:B) {
  
  
  cat("\n\n\n\n*** Iteration ", b, "\n\n\n")
  
  set.seed(b*50)
  
  
  
  # Simulate correlated covariates -----------------------------
  mysigma = matrix(1, P, P)
  
  for (i in 1:P) {
    for (j in 1:P) {
      mysigma[i, j] = 0.6^abs(i - j) + ifelse(i == j, 0, 0.1)
    }
  }
  
  
  X = MASS::mvrnorm(N, mu = rep(0, P), mysigma)
  
  Pscore <- pnorm(-0.4 + 0.3*X[, 1] + 0.2*X[, 2])
  z <- rbinom(N, 1, Pscore)
  
  summary(Pscore)
  cor(X)
  
  # Simulate POs and Y
  mu = 3 + X[, 1] + 0.8*sin(X[, 2]) + 0.7*X[, 3]*X[, 4] - X[, 5]
  ITE = 2 + 0.8*X[, 1] - 0.3*X[, 2]^2
  
  Y = mu + ITE*z + rnorm(N, 0, 1)
  
  summary(Y)
  summary(mu); summary(ITE)
  
  
  
  # Estimate Pscore with 1-hidden layer NN
  PS_nn <- nnet(x = X, y = z, size = 10, maxit = 2000, 
                decay = 0.01, trace=FALSE, abstol = 1.0e-8) 
  PS_est = PS_nn$fitted.values

  

  # Run BCF -------------------------
  mybcf <-
    SparseBCF(y = Y, z = z, x_control = X, 
              pihat = PS_est, 
              OOB = F, 
              sparse = F,
              update_interval = 5000,
              nburn = 8000, nsim = 4000)  

  
  
  # Run BCF with 5 relevant covariates only (Oracle BCF) -------------------
  my5bcf <-
    SparseBCF(y = Y, z = z, x_control = X[, 1:5], 
              pihat = PS_est, 
              OOB = F, 
              sparse = F,
              update_interval = 5000,
              nburn = 8000, nsim = 4000)
  
  
  
  # Run Sparse BCF with all 10 covariates -------------------------
  mysbcf <-
    SparseBCF(y = Y, z = z, x_control = X, 
              pihat = PS_est, 
              OOB = F, 
              sparse = T,
              update_interval = 5000,
              nburn = 8000, nsim = 4000)
  
  
  # Compute metrics
  Tau_BCF <- colMeans(mybcf$tau)
  Tau_5BCF <- colMeans(my5bcf$tau)
  Tau_SBCF <- colMeans(mysbcf$tau)

  Bias_BCF[b] = bias(Tau_BCF, ITE)
  Bias_5BCF[b] = bias(Tau_5BCF, ITE)
  Bias_SBCF[b] = bias(Tau_SBCF, ITE)
  
  PEHE_BCF[b] = PEHE(Tau_BCF, ITE)
  PEHE_5BCF[b] = PEHE(Tau_5BCF, ITE)
  PEHE_SBCF[b] = PEHE(Tau_SBCF, ITE)
  
  Cover_BCF[b] = coverage_95(mybcf$tau, ITE)
  Cover_5BCF[b] = coverage_95(my5bcf$tau, ITE)
  Cover_SBCF[b] =   coverage_95(mysbcf$tau, ITE)
  
  Split_Mu[b, ] = colMeans(mysbcf$varprb_mu)
  Split_Tau[b, ] = colMeans(mysbcf$varprb_tau)
  
}


Bias_Final = data.frame(BCF = c(mean(Bias_BCF), MC_se(Bias_BCF, B)),
                        BCF_5 = c(mean(Bias_5BCF), MC_se(Bias_5BCF, B)),
                        SparseBCF = c(mean(Bias_SBCF), MC_se(Bias_SBCF, B))
)

PEHE_Final = data.frame(BCF = c(mean(PEHE_BCF), MC_se(PEHE_BCF, B)),
                        BCF_5 = c(mean(PEHE_5BCF), MC_se(PEHE_5BCF, B)),
                        SparseBCF = c(mean(PEHE_SBCF), MC_se(PEHE_SBCF, B))
)

Cover_Final = data.frame(BCF = c(mean(Cover_BCF), MC_se(Cover_BCF, B)),
                         BCF_5 = c(mean(Cover_5BCF), MC_se(Cover_5BCF, B)),
                         SparseBCF = c(mean(Cover_SBCF), MC_se(Cover_SBCF, B))
)


Split_Proba = data.frame(Mu = colMeans(Split_Mu),
                         Tau = c(colMeans(Split_Tau), NA))


Bias_Final
PEHE_Final
Cover_Final

# Saving Results
write.csv(Bias_Final,
          paste0("YOUR DIRECTORY/Results"))

write.csv(PEHE_Final,
          paste0("YOUR DIRECTORY/Results"))

write.csv(Cover_Final,
          paste0("YOUR DIRECTORY/Results"))

write.csv(Split_Proba,
          paste0("YOUR DIRECTORY/Results"))
