###################################################
# Some other methods to compare on Simulated Data #
###################################################

rm(list = ls())


### LIBRARIES
library(tidyverse)
library(SparseBCF) # OOB Bayesian Causal Forests
library(BART) # Main package including all the version of BART
library(grf)
library(rlearner)
library(future)
availableCores() # 8 processes in total
plan(multisession)  # use future() to assign and value() function to block subsequent evaluations


### EVALUATION FUNCTIONS
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
  X[, 1:10] = qnorm(unif[, 1:10])
  X[, 11:P] = qbinom(unif[, 11:P], 1, 0.3)
  
  return(X)
  
}


### INFO
N = 1000
P = 25
B = 1000


# Store simulations true and estimated quantities
Train_PEHE_SOLS = c(NA); Test_PEHE_SOLS = c(NA)
Train_PEHE_TOLS = c(NA); Test_PEHE_TOLS = c(NA)
Train_PEHE_RLASSO = c(NA); Test_PEHE_RLASSO = c(NA)
Train_PEHE_CRF = c(NA); Test_PEHE_CRF = c(NA)
Train_PEHE_SBART = c(NA); Test_PEHE_SBART = c(NA)
Train_PEHE_TBART = c(NA); Test_PEHE_TBART = c(NA)
Train_PEHE_SDART = c(NA); Test_PEHE_SDART = c(NA)
Train_PEHE_TDART = c(NA); Test_PEHE_TDART = c(NA)
Train_PEHE_BCF = c(NA); Test_PEHE_BCF = c(NA)
Train_PEHE_SparseBCF = c(NA); Test_PEHE_SparseBCF = c(NA)


# Variable Selection
SplitProba_mu = matrix(NA, B, (P+1))  # Augmented with Propensity Score
SplitProba_tau = matrix(NA, B, P)


system.time(
  
  
  for (i in 1:B) {
    
    
    gc()
    
    cat("\n\n\n\n-------- Iteration", i, "--------\n\n\n\n")
    
    
    if(i<=500){set.seed(100 + i*5)}
    if(i>500){set.seed(5000 + i*5)}
    
    
    
    # Simulate covariates
    X <- get_features(N, P)
    
    
    # Simulate Pscore and Z
    und_lin = -0.5 + 0.2*X[, 1] + 0.1*X[, 2] + 0.4*X[, 21] + runif(N)/10
    pscore <- pnorm( und_lin )
    
    Z = rbinom(N, 1, pscore)
    
    
    
    ############## RESPONSE SURFACES
    mu = 3 + 1.5*sin(pi*X[, 1]) + 0.5*(X[, 2] - 0.5)^2 + 1.5*(2-abs(X[, 3])) + 1.5*X[, 4]*(X[, 21] + 1)
    ITE = 0.1 +  1*abs(X[, 1] - 1)*(X[, 21] + 2)
    
    sigma = sd(mu)/2
    Y = mu + ITE*Z + rnorm(N, 0, sigma)
    
    
    
    # PS estimation (we use Probit BART in this case as classifier, just to show different models are suitable)
    myBART_touse <- pbart(X, Z, sparse = FALSE, nskip = 2000, 
                          ndpost = 4000, printevery = 6000)
    
    pred_BART_touse <- myBART_touse$prob.train.mean
    sprintf("%.5f", BScore(pred_BART_touse, Z))
    
    
    # Train-Test Split
    mysplit <- c(rep(1, ceiling(0.7*N)), 
                 rep(2, floor(0.3*N)))
    
    smp_split <- sample(mysplit, replace = FALSE)  # random permutation
    
    y_train <- Y[smp_split == 1]
    y_test <- Y[smp_split == 2]
    
    X_train <- X[smp_split == 1, ]
    X_test <- X[smp_split == 2, ]
    
    z_train <- Z[smp_split == 1]
    z_test <- Z[smp_split == 2]
    
    PS_train <- pred_BART_touse[smp_split == 1]
    PS_test <- pred_BART_touse[smp_split == 2]
    
    Train_ITE <- ITE[smp_split == 1]
    Test_ITE <- ITE[smp_split == 2]
    
    
    # PS-augmented covariates
    train_augmX <- cbind(X_train, PS_train)
    test_augmX <- cbind(X_test, PS_test)
    
    colnames(train_augmX)[colnames(train_augmX) == "PS_train"] <- ""
    colnames(test_augmX)[colnames(test_augmX) == "PS_test"] <- ""
    
    
    
    
    ###### MODELS ESTIMATION  ------------------------------------------------
    
    ######################### S-Linear Regression
    # Train
    SOLS <- lm(cbind.data.frame(y_train, z_train, train_augmX))
    
    Y0_train <- predict(SOLS, newdata = cbind.data.frame( train_augmX, z_train=0) , se.fit = T)
    Y1_train <- predict(SOLS, newdata = cbind.data.frame( train_augmX, z_train=1) , se.fit = T)
    
    Train_PEHE_SOLS[i] = PEHE(Train_ITE, Y1_train$fit - Y0_train$fit)
    
    
    # Test
    Y0_test <- predict(SOLS, newdata = cbind.data.frame( test_augmX, z_train=0) , se.fit = T)
    Y1_test <- predict(SOLS, newdata = cbind.data.frame( test_augmX, z_train=1) , se.fit = T)
    
    Test_PEHE_SOLS[i] = PEHE(Test_ITE, Y1_test$fit - Y0_test$fit)
    
    
    
    ######################### T-Linear Regression
    # Train
    
    TOLS1 <- lm(cbind.data.frame(y_train, z_train, train_augmX)[z_train == 1, ])
    TOLS0 <- lm(cbind.data.frame(y_train, z_train, train_augmX)[z_train == 0, ])
    
    Y0_train <- predict(TOLS0, newdata = cbind.data.frame( train_augmX, z_train==0) , se.fit = T)
    Y1_train <- predict(TOLS1, newdata = cbind.data.frame( train_augmX, z_train==1) , se.fit = T)
    
    Train_PEHE_TOLS[i] = PEHE(Train_ITE, Y1_train$fit - Y0_train$fit)
    
    
    # Test
    Y0_test <- predict(TOLS0, newdata = cbind.data.frame( test_augmX, z_train=0) , se.fit = T)
    Y1_test <- predict(TOLS1, newdata = cbind.data.frame( test_augmX, z_train=1) , se.fit = T)
    
    Test_PEHE_TOLS[i] = PEHE(Test_ITE, Y1_test$fit - Y0_test$fit)
    
    
    
    
    ######################## R-Lasso Regression
    # No estimated PS as 
    RLASSO <- rlasso(x = train_augmX[, -ncol(train_augmX)], w = z_train, y = y_train)
    
    Train_RL_tau = predict(RLASSO, train_augmX[, -ncol(train_augmX)])
    Test_RL_tau = predict(RLASSO, test_augmX[, -ncol(train_augmX)])
    
    Train_PEHE_RLASSO[i] = PEHE(Train_ITE, Train_RL_tau)
    Test_PEHE_RLASSO[i] = PEHE(Test_ITE, Test_RL_tau)
    
    rm(RLASSO)
    
    
    
    ######################### Causal RF (Athey, Wagner 2016)
    # Train
    CRF <- causal_forest(train_augmX, y_train, z_train)
    
    # Test
    Train_PEHE_CRF[i] = PEHE(Train_ITE, predict(CRF)[, "predictions"])
    Test_PEHE_CRF[i] = PEHE(Test_ITE, predict(CRF, newdata =  test_augmX)[, "predictions"])
    
    
    # Remove garbage
    rm(CRF)
    
    
    
    ######################### S-BART
    #### Train
    myBART <- wbart(x.train = cbind(train_augmX, z_train), y.train = y_train, 
                    x.test = cbind(test_augmX, z_test), nskip = 2000, ndpost = 4000, printevery = 6000)
    
    XZ0_train <- cbind(train_augmX, z_train)
    XZ0_train[, "z_train"] <- ifelse(XZ0_train[, "z_train"] == 1, 0, 1)
    
    
    Y0_train <- predict(myBART, newdata = XZ0_train)
    
    
    # Non Sparse
    All_obs <- cbind(Y1 = myBART$yhat.train.mean, 
                     train_augmX, 
                     z_train)
    All_count <- cbind(Y0 = colMeans(Y0_train), 
                       XZ0_train)
    
    All_Trt <- All_obs
    All_Trt[which(All_Trt[, "z_train"] == 0), ] <- All_count[which(All_count[, "z_train"] == 1), ]
    
    All_Ctrl <- All_count
    All_Ctrl[which(All_Ctrl[, "z_train"] == 1), ] <- All_obs[which(All_obs[, "z_train"] == 0), ]
    
    # Store estimates
    Train_PEHE_SBART[i] = PEHE(Train_ITE, All_Trt[, "Y1"] - All_Ctrl[, "Y0"])
    
    
    
    #### Test
    XZ0_test <- cbind(test_augmX, z_test)
    XZ0_test[, "z_test"] <- ifelse(XZ0_test[, "z_test"] == 1, 0, 1)
    
    Y0_test <- predict(myBART, newdata = XZ0_test)
    
    
    # Non Sparse
    All_obs <- cbind(Y1 = myBART$yhat.test.mean, 
                     test_augmX, 
                     z_test)
    All_count <- cbind(Y0 = colMeans(Y0_test), 
                       XZ0_test)
    
    All_Trt <- All_obs
    All_Trt[which(All_Trt[, "z_test"] == 0), ] <- All_count[which(All_count[, "z_test"] == 1), ]
    
    All_Ctrl <- All_count
    All_Ctrl[which(All_Ctrl[, "z_test"] == 1), ] <- All_obs[which(All_obs[, "z_test"] == 0), ]
    
    # Store estimates
    Test_PEHE_SBART[i] = PEHE(Test_ITE, All_Trt[, "Y1"] - All_Ctrl[, "Y0"])
    
    
    
    # Free up space
    rm(myBART, Y0_train, Y0_test)
    
    
    
    ###########
    # Check that there is no single-valued binary variable and in case remove it for T-BART (as it would not split on it)
    check_train1 = which(apply(train_augmX[z_train == 1, ], 2, function(x) max(x) == min(x)))
    check_train0 = which(apply(train_augmX[z_train == 0, ], 2, function(x) max(x) == min(x)))
    
    check_test1 = which(apply(test_augmX[z_test == 1, ], 2, function(x) max(x) == min(x)))
    check_test0 = which(apply(test_augmX[z_test == 0, ], 2, function(x) max(x) == min(x)))
    
    check_all = c(check_train1, check_train0, check_test1, check_test0)
    
    if (length(check_all) > 0) {
      
      new_train_X = train_augmX[, -check_all]
      new_test_X = test_augmX[, -check_all]
      
    } else {
      
      new_train_X = train_augmX
      new_test_X = test_augmX
      
    }
    
    
    ######################### T-BART
    #### Train
    myBART1_shell <- future({
      wbart(x.train = new_train_X[z_train == 1, ], y.train = y_train[z_train == 1], x.test = new_test_X[z_test == 1, ],
            nskip = 2000, ndpost = 4000, printevery = 6000)
    })
    
    myBART0_shell <- future({
      wbart(x.train = new_train_X[z_train == 0, ], y.train = y_train[z_train == 0], x.test = new_test_X[z_test == 0, ],
            nskip = 2000, ndpost = 4000, printevery = 6000)
    })  
    
    myBART1 <- value(myBART1_shell); myBART0 <- value(myBART0_shell)
    rm(myBART1_shell, myBART0_shell)
    
    
    # (predict counterfactual)
    Y1_0_shell <- future({
      predict(myBART1, newdata = new_train_X[z_train == 0, ])
    })
    
    Y0_1 <- predict(myBART0, newdata = new_train_X[z_train == 1, ])
    
    Y1_0 <- value(Y1_0_shell)
    rm(Y1_0_shell)
    
    Y1_train <- y_train
    Y1_train[z_train == 1] <- myBART1$yhat.train.mean
    Y1_train[z_train == 0] <- colMeans(Y1_0)
    
    Y0_train <- y_train
    Y0_train[z_train == 0] <- myBART0$yhat.train.mean
    Y0_train[z_train == 1] <- colMeans(Y0_1)
    
    # Store ITE estimates
    Train_PEHE_TBART[i] = PEHE(Train_ITE, Y1_train - Y0_train)
    
    
    
    ##### Test
    
    # (predict counterfactual)
    Y1_0_shell <- future({
      predict(myBART1, newdata = new_test_X[z_test == 0, ])
    })
    
    Y0_1 <- predict(myBART0, newdata = new_test_X[z_test == 1, ])
    
    Y1_0 <- value(Y1_0_shell)
    rm(Y1_0_shell)
    
    Y1_test <- y_test
    Y1_test[z_test == 1] <- myBART1$yhat.test.mean
    Y1_test[z_test == 0] <- colMeans(Y1_0)
    
    Y0_test <- y_test
    Y0_test[z_test == 0] <- myBART0$yhat.test.mean
    Y0_test[z_test == 1] <- colMeans(Y0_1)
    
    # Store ITE estimates
    Test_PEHE_TBART[i] = PEHE(Test_ITE, Y1_test - Y0_test)
    
    
    # Remove garbage
    rm(myBART1, myBART0, Y0_1, Y1_0, new_train_X, new_test_X)
    
    
    
    
    ######################### S-DART
    #### Train
    myDART <- wbart(x.train = cbind(train_augmX, z_train), y.train = y_train, sparse = TRUE,
                    x.test = cbind(test_augmX, z_test), nskip = 2000, ndpost = 4000, printevery = 6000)
    
    XZ0_train <- cbind(train_augmX, z_train)
    XZ0_train[, "z_train"] <- ifelse(XZ0_train[, "z_train"] == 1, 0, 1)
    
    
    Y0_train <- predict(myDART, newdata = XZ0_train)
    
    
    # Non Sparse
    All_obs <- cbind(Y1 = myDART$yhat.train.mean, 
                     train_augmX, 
                     z_train)
    All_count <- cbind(Y0 = colMeans(Y0_train), 
                       XZ0_train)
    
    All_Trt <- All_obs
    All_Trt[which(All_Trt[, "z_train"] == 0), ] <- All_count[which(All_count[, "z_train"] == 1), ]
    
    All_Ctrl <- All_count
    All_Ctrl[which(All_Ctrl[, "z_train"] == 1), ] <- All_obs[which(All_obs[, "z_train"] == 0), ]
    
    # Store estimates
    Train_PEHE_SDART[i] = PEHE(Train_ITE, All_Trt[, "Y1"] - All_Ctrl[, "Y0"])
    
    
    
    #### Test
    XZ0_test <- cbind(test_augmX, z_test)
    XZ0_test[, "z_test"] <- ifelse(XZ0_test[, "z_test"] == 1, 0, 1)
    
    Y0_test <- predict(myDART, newdata = XZ0_test)
    
    
    # Non Sparse
    All_obs <- cbind(Y1 = myDART$yhat.test.mean, 
                     test_augmX, 
                     z_test)
    All_count <- cbind(Y0 = colMeans(Y0_test), 
                       XZ0_test)
    
    All_Trt <- All_obs
    All_Trt[which(All_Trt[, "z_test"] == 0), ] <- All_count[which(All_count[, "z_test"] == 1), ]
    
    All_Ctrl <- All_count
    All_Ctrl[which(All_Ctrl[, "z_test"] == 1), ] <- All_obs[which(All_obs[, "z_test"] == 0), ]
    
    # Store estimates
    Test_PEHE_SDART[i] = PEHE(Test_ITE, All_Trt[, "Y1"] - All_Ctrl[, "Y0"])
    
    
    
    # Free up space
    rm(myDART, Y0_train, Y0_test)
    
    
    
    
    
    ###########
    # Check that there is no single-valued binary variable and in case remove it for T-DART (as it would not split on it)
    check_train1 = which(apply(train_augmX[z_train == 1, ], 2, function(x) max(x) == min(x)))
    check_train0 = which(apply(train_augmX[z_train == 0, ], 2, function(x) max(x) == min(x)))
    
    check_test1 = which(apply(test_augmX[z_test == 1, ], 2, function(x) max(x) == min(x)))
    check_test0 = which(apply(test_augmX[z_test == 0, ], 2, function(x) max(x) == min(x)))
    
    check_all = c(check_train1, check_train0, check_test1, check_test0)
    
    if (length(check_all) > 0) {
      
      new_train_X = train_augmX[, -check_all]
      new_test_X = test_augmX[, -check_all]
      
    } else {
      
      new_train_X = train_augmX
      new_test_X = test_augmX
      
    }
    
    
    ######################### T-DART
    #### Train
    myDART1_shell <- future({
      wbart(x.train = new_train_X[z_train == 1, ], y.train = y_train[z_train == 1], x.test = new_test_X[z_test == 1, ],
            nskip = 2000, ndpost = 4000, printevery = 6000, sparse = TRUE)
    })
    
    myDART0_shell <- future({
      wbart(x.train = new_train_X[z_train == 0, ], y.train = y_train[z_train == 0], x.test = new_test_X[z_test == 0, ],
            nskip = 2000, ndpost = 4000, printevery = 6000, sparse = TRUE)
    })  
    
    myDART1 <- value(myDART1_shell); myDART0 <- value(myDART0_shell)
    rm(myDART1_shell, myDART0_shell)
    
    
    # (predict counterfactual)
    Y1_0_shell <- future({
      predict(myDART1, newdata = new_train_X[z_train == 0, ])
    })
    
    Y0_1 <- predict(myDART0, newdata = new_train_X[z_train == 1, ])
    
    Y1_0 <- value(Y1_0_shell)
    rm(Y1_0_shell)
    
    Y1_train <- y_train
    Y1_train[z_train == 1] <- myDART1$yhat.train.mean
    Y1_train[z_train == 0] <- colMeans(Y1_0)
    
    Y0_train <- y_train
    Y0_train[z_train == 0] <- myDART0$yhat.train.mean
    Y0_train[z_train == 1] <- colMeans(Y0_1)
    
    # Store ITE estimates
    Train_PEHE_TDART[i] = PEHE(Train_ITE, Y1_train - Y0_train)
    
    
    
    ##### Test
    
    # (predict counterfactual)
    Y1_0_shell <- future({
      predict(myDART1, newdata = new_test_X[z_test == 0, ])
    })
    
    Y0_1 <- predict(myDART0, newdata = new_test_X[z_test == 1, ])
    
    Y1_0 <- value(Y1_0_shell)
    rm(Y1_0_shell)
    
    Y1_test <- y_test
    Y1_test[z_test == 1] <- myDART1$yhat.test.mean
    Y1_test[z_test == 0] <- colMeans(Y1_0)
    
    Y0_test <- y_test
    Y0_test[z_test == 0] <- myDART0$yhat.test.mean
    Y0_test[z_test == 1] <- colMeans(Y0_1)
    
    # Store ITE estimates
    Test_PEHE_TDART[i] = PEHE(Test_ITE, Y1_test - Y0_test)
    
    
    # Remove garbage
    rm(myDART1, myDART0, Y0_1, Y1_0, new_train_X, new_test_X)
    
    
    
    
    
    ######################### Normal BCF
    #### Train
    mybcf <-
      SparseBCF(y = y_train, 
                z = z_train, 
                x_control = X_train, 
                pihat = PS_train, 
                OOB = T, 
                sparse = F,
                x_pred_mu = X_test, 
                pi_pred = PS_test, 
                x_pred_tau = X_test,
                update_interval = 5000, 
                nburn = 3000, nsim = 7000)
    
    
    Train_PEHE_BCF[i] = PEHE(Train_ITE, apply(mybcf$tau, 2, mean))
    Test_PEHE_BCF[i] = PEHE(Test_ITE, apply(mybcf$tau_pred, 2, mean))
    
    # Remove garbage
    rm(mybcf)
    
  
    
    ######################### Sparse BCF
    #### Train
    mysbcf <-
      SparseBCF(y = y_train, 
                z = z_train, 
                x_control = X_train, 
                pihat = PS_train, 
                OOB = T, 
                sparse = T,
                x_pred_mu = X_test, 
                pi_pred = PS_test, 
                x_pred_tau = X_test,
                update_interval = 5000,
                nburn = 3000, nsim = 7000)
    
    
    Train_PEHE_SparseBCF[i] = PEHE(Train_ITE, apply(mysbcf$tau, 2, mean))
    Test_PEHE_SparseBCF[i] = PEHE(Test_ITE, apply(mysbcf$tau_pred, 2, mean))
    
    
    # Save posterior splitting probabilities
    SplitProba_mu[i, ] = colMeans(mysbcf$varprb_mu)
    SplitProba_tau[i, ] = colMeans(mysbcf$varprb_tau)
    
   
    # Remove garbage
    rm(mysbcf)
 
    
  }
  
)




# Save Results --------------------------------------------------
# PEHE
PEHE_Final = data.frame(SOLS = c(mean(Train_PEHE_SOLS), MC_se(Train_PEHE_SOLS, B),
                                 mean(Test_PEHE_SOLS), MC_se(Test_PEHE_SOLS, B)),
                        TOLS = c(mean(Train_PEHE_TOLS), MC_se(Train_PEHE_TOLS, B),
                                 mean(Test_PEHE_TOLS), MC_se(Test_PEHE_TOLS, B)),
                        RLASSO = c(mean(Train_PEHE_RLASSO), MC_se(Train_PEHE_RLASSO, B),
                                   mean(Test_PEHE_RLASSO), MC_se(Test_PEHE_RLASSO, B)),
                        CRF = c(mean(Train_PEHE_CRF), MC_se(Train_PEHE_CRF, B),
                                mean(Test_PEHE_CRF), MC_se(Test_PEHE_CRF, B)),
                        SBART = c(mean(Train_PEHE_SBART), MC_se(Train_PEHE_SBART, B),
                                  mean(Test_PEHE_SBART), MC_se(Test_PEHE_SBART, B)),
                        TBART = c(mean(Train_PEHE_TBART), MC_se(Train_PEHE_TBART, B),
                                  mean(Test_PEHE_TBART), MC_se(Test_PEHE_TBART, B)),
                        SDART = c(mean(Train_PEHE_SDART), MC_se(Train_PEHE_SDART, B),
                                  mean(Test_PEHE_SDART), MC_se(Test_PEHE_SDART, B)),
                        TDART = c(mean(Train_PEHE_TDART), MC_se(Train_PEHE_TDART, B),
                                  mean(Test_PEHE_TDART), MC_se(Test_PEHE_TDART, B)),
                        BCF = c(mean(Train_PEHE_BCF), MC_se(Train_PEHE_BCF, B),
                                mean(Test_PEHE_BCF), MC_se(Test_PEHE_BCF, B)),
                        SparseBCF = c(mean(Train_PEHE_SparseBCF), MC_se(Train_PEHE_SparseBCF, B),
                                      mean(Test_PEHE_SparseBCF), MC_se(Test_PEHE_SparseBCF, B))
)
rownames(PEHE_Final) = c("Train", "SE_Train", "Test", "SE_Test")



# All iterations results
PEHE_Single = data.frame(SOLS_Train = Train_PEHE_SOLS, SOLS_Test = Test_PEHE_SOLS,
                         TOLS_Train = Train_PEHE_TOLS, TOLS_Test = Test_PEHE_TOLS,
                         RLASSO_Train = Train_PEHE_RLASSO, RLASSO_Test = Test_PEHE_RLASSO,
                         CRF_Train = Train_PEHE_CRF, CRF_Test = Test_PEHE_CRF,
                         SBART_Train = Train_PEHE_SBART, SBART_Test = Test_PEHE_SBART,
                         TBART_Train = Train_PEHE_TBART, TBART_Test = Test_PEHE_TBART,
                         SDART_Train = Train_PEHE_SDART, SDART_Test = Test_PEHE_SDART,
                         TDART_Train = Train_PEHE_TDART, TDART_Test = Test_PEHE_TDART,
                         BCF_Train = Train_PEHE_BCF, BCF_Test = Test_PEHE_BCF,
                         SparseBCF_Train = Train_PEHE_SparseBCF, SparseBCF_Test = Test_PEHE_SparseBCF)



# Posterior Splitting Probabilities
barplot(colMeans(SplitProba_mu))
barplot(colMeans(SplitProba_tau))

Split_Final = data.frame(Mu = colMeans(SplitProba_mu),
                         Tau = c(colMeans(SplitProba_tau), NA)
)



# Save PEHE
curr_dir <- dirname(rstudioapi::getSourceEditorContext()$path)

write.csv(PEHE_Final,
          paste0(curr_dir, "/Results/PEHE_Final_P25.csv"))


write.csv(PEHE_Single,
          paste0(curr_dir, "/Results/PEHE_Single_P25.csv"))


write.csv(Split_Final,
          paste0(curr_dir, "/Results/Split_Final_P25.csv"))

