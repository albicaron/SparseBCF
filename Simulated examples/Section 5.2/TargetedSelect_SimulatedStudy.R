##########################################
# Simulated example of nudging PS splits #
##########################################

rm(list = ls())


### LIBRARIES
library(tidyverse)
library(SparseBCF) # OOB Sp. Bayesian Causal Forests
library(BART) # Main package including all the version of BART
library(nnet)


### FUNCTIONS
# Evaluatiion
bias <- function(x,y) mean(x-y)
PEHE <- function(x,y) sqrt(mean((x-y)^2))
coverage_95 <- function(ITE_EST, ITE) {
  quan = apply(ITE_EST, 2, function(x) quantile(x, c(0.025, 0.975)))
  cove = sum((quan[1,] < ITE) & (ITE < quan[2,]))/length(ITE)
  
  return(cove)
}
MC_se <- function(x, B) qt(0.975, B-1)*sd(x)/sqrt(B)


# Correlated covariates
get_features <- function(N, P) {
  
  # Generate correlated uniforms from a Gaussian Copula
  mysigma = matrix(1, P, P)
  
  for (i in 1:P) {
    for (j in 1:P) {
      mysigma[i, j] = 0.5^abs(i - j) + ifelse(i == j, 0, 0.1)
    }
  }
  
  mycop = MASS::mvrnorm(N, rep(0, P), Sigma = mysigma)
  unif = pnorm(mycop)
  
  
  # Transform in continuous and binary covariates
  X = matrix(NA, N, P)
  X[, 1:5] = qnorm(unif[, 1:5])
  X[, 6:P] = qbinom(unif[, 6:P], 1, 0.3)
  
  return(X)
}


### INFO
N = 1000
P = 10  # We try P = 10 and P = 20 to see how the problem degenerate
B = 250

# Simulate correlated covariates
X <- get_features(N, P)


# Store Metrics
Bias_BCF = c(NA); Bias_SBCF = c(NA); Bias_PSBCFA = c(NA); Bias_PSBCFB = c(NA)
PEHE_BCF = c(NA); PEHE_SBCF = c(NA); PEHE_PSBCFA = c(NA); PEHE_PSBCFB = c(NA)
Cover_BCF = c(NA); Cover_SBCF = c(NA); Cover_PSBCFA = c(NA); Cover_PSBCFB = c(NA)

Bias_noPS = c(NA) 
PEHE_noPS = c(NA)
Cover_noPS = c(NA)


Split_Mu_SBCF = matrix(NA, B, P+1); Split_Tau_SBCF = matrix(NA, B, P)
Split_Mu_PSBCFA = matrix(NA, B, P+1); Split_Tau_PSBCFA = matrix(NA, B, P)
Split_Mu_PSBCFB = matrix(NA, B, P+1); Split_Tau_PSBCFB = matrix(NA, B, P)
Split_Mu_noPS = matrix(NA, B, P); Split_Tau_noPS = matrix(NA, B, P)


# Simulation Study start
system.time(
  
  for (b in 1:B) {
    
    gc()
    cat("\n\n\n\n*** Iteration ", b, "\n\n\n")
    
    set.seed(b*5)
    
    
    # Simulate correlated covariates ---------------------------------
    X <- get_features(N, P)
    
    # Generate quantities, with strong targeted selection
    mu_base = 2 + 0.5*sin(pi*X[, 1]) - 0.25*X[, 2]^2 + 0.75*X[, 3]*X[, 9]  
    mu = 5*mu_base
    # plot(X[, 1],  mu)    # Uncomment to check quantities
    
    ITE = 1 + 2*abs(X[, 4]) + 1*X[, 10]
    # plot(X[, 4], ITE)           # Uncomment to check quantities
    # summary(mu); summary(ITE)    
    # par(mfrow = c(1, 2))
    # hist(mu); hist(ITE)
    
    Pscore = 0.9*plogis(1.2 - mu_base)
    z <- rbinom(N, 1, Pscore)
    # summary(Pscore)   # Uncomment to check quantities
    # hist(Pscore)
    # table(z)
    
    Y = mu + ITE*z + rnorm(N, 0, sd(ITE)/2)
    # summary(Y)   # Uncomment to check quantities
    # summary(mu); summary(ITE)
    # 
    
    
    # We do not estimate PS, but assume we know it fully  --------------------------
    
    
    # Default BCF  -------------------------
    mybcf <-
      SparseBCF(y = Y, z = z, x_control = X, 
                pihat = Pscore, 
                OOB = F, 
                sparse = F,
                update_interval = 5000,
                nburn = 8000, nsim = 4000)
    
    
    # Sparse BCF -------------------------
    mysbcf <-
      SparseBCF(y = Y, z = z, x_control = X, 
                pihat = Pscore, 
                OOB = F, 
                sparse = T,
                update_interval = 5000,
                nburn = 8000, nsim = 4000)
    
    
    # NO-PS Sparse BCF -------------------------
    # To implement this we circumvent the SparseBCF plugging the last covariate as PS 
    sbcf_noPS <-
      SparseBCF(y = Y, z = z, x_control = X[, -ncol(X)], x_moderate = X,
                pihat = X[, ncol(X)], 
                OOB = F, 
                sparse = T,
                update_interval = 5000,
                nburn = 8000, nsim = 4000)
    
    
    # Informative Prior BCF with splits on the PS -------------------------
    # Only \mu has an informed prior in this case;
    # While we remain agnostic about \tau
    myPSbcfA <-
      SparseBCF(y = Y, z = z, x_control = X, 
                pihat = Pscore, 
                OOB = F, 
                sparse = T, 
                inform_mu = T, weights_mu = c(rep(1, P), 50),     # PS is at the last spot
                inform_tau = F,
                update_interval = 5000,
                nburn = 8000, nsim = 4000)
    
    
    myPSbcfB <-
      SparseBCF(y = Y, z = z, x_control = X, 
                pihat = Pscore, 
                OOB = F, 
                sparse = T, 
                inform_mu = T, weights_mu = c(rep(1, P), 100),     # PS is at the last spot
                inform_tau = F,
                update_interval = 5000,
                nburn = 8000, nsim = 4000)
    
    
    # Compute metrics
    Tau_BCF <- colMeans(mybcf$tau)
    Tau_SBCF <- colMeans(mysbcf$tau)
    Tau_noPS <- colMeans(sbcf_noPS$tau)
    Tau_PSBCFA <- colMeans(myPSbcfA$tau)
    Tau_PSBCFB <- colMeans(myPSbcfB$tau)
    
    Bias_BCF[b] = bias(Tau_BCF, ITE)
    Bias_noPS[b] = bias(Tau_noPS, ITE)
    Bias_SBCF[b] = bias(Tau_SBCF, ITE)
    Bias_PSBCFA[b] = bias(Tau_PSBCFA, ITE)
    Bias_PSBCFB[b] = bias(Tau_PSBCFB, ITE)
    
    PEHE_BCF[b] = PEHE(Tau_BCF, ITE)
    PEHE_SBCF[b] = PEHE(Tau_SBCF, ITE)
    PEHE_noPS[b] = PEHE(Tau_noPS, ITE)
    PEHE_PSBCFA[b] = PEHE(Tau_PSBCFA, ITE)
    PEHE_PSBCFB[b] = PEHE(Tau_PSBCFB, ITE)
    
    Cover_BCF[b] = coverage_95(mybcf$tau, ITE)
    Cover_SBCF[b] = coverage_95(mysbcf$tau, ITE)
    Cover_noPS[b] = coverage_95(sbcf_noPS$tau, ITE)
    Cover_PSBCFA[b] = coverage_95(myPSbcfA$tau, ITE)
    Cover_PSBCFB[b] = coverage_95(myPSbcfB$tau, ITE)
    
    # Posterior splitting probs
    Split_Mu_SBCF[b, ] = colMeans(mysbcf$varprb_mu)
    Split_Tau_SBCF[b, ] = colMeans(mysbcf$varprb_tau)
    
    Split_Mu_noPS[b, ] = colMeans(sbcf_noPS$varprb_mu)
    Split_Tau_noPS[b, ] = colMeans(sbcf_noPS$varprb_tau)
    
    Split_Mu_PSBCFA[b, ] = colMeans(myPSbcfA$varprb_mu)
    Split_Tau_PSBCFA[b, ] = colMeans(myPSbcfA$varprb_tau)
    
    Split_Mu_PSBCFB[b, ] = colMeans(myPSbcfB$varprb_mu)
    Split_Tau_PSBCFB[b, ] = colMeans(myPSbcfB$varprb_tau)
    
    
  }
)



# Compute averaged quantities
Bias_Final = data.frame(BCF = c(mean(Bias_BCF), MC_se(Bias_BCF, B)),
                        SparseBCF = c(mean(Bias_SBCF), MC_se(Bias_SBCF, B)),
                        SparseBCF_noPS = c(mean(Bias_noPS), MC_se(Bias_noPS, B)),
                        PS_BCF50 = c(mean(Bias_PSBCFA), MC_se(Bias_PSBCFA, B)),
                        PS_BCF100 = c(mean(Bias_PSBCFB), MC_se(Bias_PSBCFB, B))
                        
)

PEHE_Final = data.frame(BCF = c(mean(PEHE_BCF), MC_se(PEHE_BCF, B)),
                        SparseBCF = c(mean(PEHE_SBCF), MC_se(PEHE_SBCF, B)),
                        SparseBCF_noPS = c(mean(PEHE_noPS), MC_se(PEHE_noPS, B)),
                        PS_BCF50 = c(mean(PEHE_PSBCFA), MC_se(PEHE_PSBCFA, B)),
                        PS_BCF100 = c(mean(PEHE_PSBCFB), MC_se(PEHE_PSBCFB, B))
)

Cover_Final = data.frame(BCF = c(mean(Cover_BCF), MC_se(Cover_BCF, B)),
                         SparseBCF = c(mean(Cover_SBCF), MC_se(Cover_SBCF, B)),
                         SparseBCF_noPS = c(mean(Cover_noPS), MC_se(Cover_noPS, B)),
                         PS_BCF50 = c(mean(Cover_PSBCFA), MC_se(Cover_PSBCFA, B)),
                         PS_BCF100 = c(mean(Cover_PSBCFB), MC_se(Cover_PSBCFB, B))
)



Split_Proba = data.frame(Mu_SBCF = colMeans(Split_Mu_SBCF),
                         Tau_SBCF = c(colMeans(Split_Tau_SBCF), NA),
                         Mu_SBCF_noPS = c(colMeans(Split_Mu_noPS), NA),
                         Tau_SBCF_noPS = c(colMeans(Split_Tau_noPS), NA),
                         Mu_BCF50 = colMeans(Split_Mu_PSBCFA),
                         Tau_BCF50 = c(colMeans(Split_Tau_PSBCFA), NA),
                         Mu_BCF100 = colMeans(Split_Mu_PSBCFB),
                         Tau_BCF100 = c(colMeans(Split_Tau_PSBCFB), NA))


Bias_Final
PEHE_Final
Cover_Final


write.csv(Bias_Final,
          paste0("YOUR DIRECTORY/TruePSversion_ConfExam_FinalBias_B", B, "_P", P, ".csv"))

write.csv(PEHE_Final,
          paste0("YOUR DIRECTORY/TruePSversion_ConfExam_FinalPEHE_B", B, "_P", P, ".csv"))

write.csv(Cover_Final,
          paste0("YOUR DIRECTORY/TruePSversion_ConfExam_FinalCoverage_B", B, "_P", P, ".csv"))

write.csv(Split_Proba,
          paste0("YOUR DIRECTORY/TruePSversion_ConfExam_SplitProba_B", B, "_P", P, ".csv"))

