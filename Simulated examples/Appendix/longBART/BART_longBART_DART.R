#########################################
# Test: BART vs long-chain BART vs DART #
#########################################

rm(list = ls())

library(tidyverse)
library(BART) # Main package including all the version of BART

# Functions
RMSE <- function(x,y) sqrt(mean((x-y)^2))
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
  X[, 1:10] = qnorm(unif[, 1:10])
  X[, 11:P] = qbinom(unif[, 11:P], 1, 0.3)
  
  return(X)
}


### INFO
N = 500
P = 50
B = 500

# Store simulations true and estimated quantities
Train_RMSE_BART = c(NA); Test_RMSE_BART = c(NA)
Train_RMSE_BART_100K = c(NA); Test_RMSE_BART_100K = c(NA)
Train_RMSE_DART = c(NA); Test_RMSE_DART = c(NA)


# Splitting Counts
SplitCounts_BART = matrix(NA, B, P)
SplitCounts_BART_100K = matrix(NA, B, P)
SplitCounts_DART = matrix(NA, B, P)


system.time(

for (i in 1:B) {
  
  gc()
  cat("\n*** Iteration ", i, "\n")
  set.seed(i*11)
  
  # Simulate predictors and outcomes
  X <- get_features(N, P)
  mu <- 5 + 5*sin(pi*X[, 1]) + 2.5*(X[, 2] - 0.5)^2 + 1.5*abs(X[, 3]) + 2*X[, 4]*(X[, 20] + 1)
  
  Y <- rnorm(N, mean = mu, sd = 1)
  
  # Train-Test Split
  mysplit <- c(rep(1, ceiling(0.7*N)), 
               rep(2, floor(0.3*N)))
  
  smp_split <- sample(mysplit, replace = FALSE)  # random permutation
  
  y_train <- Y[smp_split == 1]
  y_test <- Y[smp_split == 2]
  
  X_train <- X[smp_split == 1, ]
  X_test <- X[smp_split == 2, ]
  
  # BART 10,000K
  BART = wbart(x.train = X_train, y.train = y_train, 
               x.test = X_test, sparse = FALSE, nskip = 4000, 
               ndpost = 2000, printevery = 6000)
  
  longBART = wbart(x.train = X_train, y.train = y_train, 
                   x.test = X_test, sparse = FALSE, 
                   nskip = 40000, ndpost = 20000, 
                   printevery = 50000)
  
  DART = wbart(x.train = X_train, y.train = y_train, 
               x.test = X_test, sparse = TRUE, nskip = 4000, 
               ndpost = 2000, printevery = 6000)

  
  # Evaluate
  Train_RMSE_BART[i] <- RMSE(y_train, BART$yhat.train.mean)
  Test_RMSE_BART[i] <- RMSE(y_test, BART$yhat.test.mean)
  SplitCounts_BART[i, ] <- colMeans(BART$varcount)
  
  Train_RMSE_BART_100K[i] <- RMSE(y_train, longBART$yhat.train.mean)
  Test_RMSE_BART_100K[i] <- RMSE(y_test, longBART$yhat.test.mean)
  SplitCounts_BART_100K[i, ] <- colMeans(longBART$varcount)
  
  Train_RMSE_DART[i] <- RMSE(y_train, DART$yhat.train.mean)
  Test_RMSE_DART[i] <- RMSE(y_test, DART$yhat.test.mean)
  SplitCounts_DART[i, ] <- colMeans(DART$varcount)
  
}
)

RMSE_final = data.frame(BART = c(mean(Train_RMSE_BART), MC_se(Train_RMSE_BART, B),
                                 mean(Test_RMSE_BART), MC_se(Test_RMSE_BART, B)),
                        longBART = c(mean(Train_RMSE_BART_100K), MC_se(Train_RMSE_BART_100K, B),
                                     mean(Test_RMSE_BART_100K), MC_se(Test_RMSE_BART_100K, B)),
                        DART = c(mean(Train_RMSE_DART), MC_se(Train_RMSE_DART, B),
                                 mean(Test_RMSE_DART), MC_se(Test_RMSE_DART, B)))
rownames(RMSE_final) = c("Train", "SE_Train", "Test", "SE_Test")


SplitCounts_final = data.frame(BART = colMeans(SplitCounts_BART),
                               BART_std_err = apply(SplitCounts_BART, 2, function(x) MC_se(x, B)),
                               longBART = colMeans(SplitCounts_BART_100K),
                               longBART_std_err = apply(SplitCounts_BART_100K, 2, function(x) MC_se(x, B)),
                               DART = colMeans(SplitCounts_DART),
                               DART_std_err = apply(SplitCounts_DART, 2, function(x) MC_se(x, B)))

RMSE_final
SplitCounts_final


# Save 
curr_dir <- dirname(rstudioapi::getSourceEditorContext()$path)

write.csv(RMSE_final,
          paste0(curr_dir, "/Results/RMSE_Final_P20.csv"))

write.csv(SplitCounts_final,
          paste0(curr_dir, "/Results/SplitCounts_Final_P20.csv"))


