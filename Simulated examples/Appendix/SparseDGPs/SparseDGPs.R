###################################################
# Some other methods to compare on Simulated Data #
###################################################

rm(list = ls())


### LIBRARIES
library(tidyverse)
library(SparseBCF) # OOB Bayesian Causal Forests
library(BART) # Main package including all the version of BART

### EVALUATION FUNCTIONS
curr_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
source(paste0(curr_dir, "/utils.R"))

dgp_list = list(NOsp_dgp, sp_PS_dgp, sp_Mu_dgp, sp_Tau_dgp)

### INFO
N = 500
P = 5
B = 200


# Store simulations true and estimated quantities
Train_PEHE_BCF = matrix(NA, B, length(dgp_list)); Test_PEHE_BCF = matrix(NA, B, length(dgp_list))
Train_PEHE_SparseBCF = matrix(NA, B, length(dgp_list)); Test_PEHE_SparseBCF = matrix(NA, B, length(dgp_list))

for (j in 1:length(dgp_list)) {
  for (i in 1:B) {
    
    gc()
    
    cat("\n\n\n\n-------- Iteration", i, "--------\n\n\n\n")
    
    
    if(i<=500){set.seed(100 + i*5)}
    if(i>500){set.seed(5000 + i*5)}
    
    
    # Simulate DGP
    X <- get_features(N, P)
    mylist = dgp_list[[j]](N, X)
    
    Y = mylist$Y
    Z = mylist$Z
    ITE = mylist$ITE
    
    rm(mylist)
    
    # PS estimation (we use Probit BART in this case as classifier, just to show different models are suitable)
    myBART_touse <- pbart(X, Z, sparse = FALSE, nskip = 2000, 
                          ndpost = 4000, printevery = 6000)
    
    pred_BART_touse <- myBART_touse$prob.train.mean
    # sprintf("%.5f", BScore(pred_BART_touse, Z))
    
    
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
    
    
    Train_PEHE_BCF[i, j] = PEHE(Train_ITE, apply(mybcf$tau, 2, mean))
    Test_PEHE_BCF[i, j] = PEHE(Test_ITE, apply(mybcf$tau_pred, 2, mean))
    
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
    
    
    Train_PEHE_SparseBCF[i, j] = PEHE(Train_ITE, apply(mysbcf$tau, 2, mean))
    Test_PEHE_SparseBCF[i, j] = PEHE(Test_ITE, apply(mysbcf$tau_pred, 2, mean))
    
    # Remove garbage
    rm(mysbcf)
    
  }
}



# Save Results --------------------------------------------------
# PEHE
colMeans(Train_PEHE_BCF); colMeans(Test_PEHE_BCF)
colMeans(Train_PEHE_SparseBCF); colMeans(Test_PEHE_SparseBCF)

df_tr = data.frame(PEHE_BCF = colMeans(Train_PEHE_BCF), 
                   SE_BCF = apply(Train_PEHE_BCF, 2, function(x) MC_se(x, B)),
                   PEHE_SHBCF = colMeans(Train_PEHE_SparseBCF), 
                   SE_SHBCF = apply(Train_PEHE_SparseBCF, 2, function(x) MC_se(x, B)),
                   DGP = c("Not Sparse", "Sparse PS", "Sparse Mu", "Sparse Tau"))

df_te = data.frame(PEHE_BCF = colMeans(Test_PEHE_BCF), 
                   SE_BCF = apply(Test_PEHE_BCF, 2, function(x) MC_se(x, B)),
                   PEHE_SHBCF = colMeans(Test_PEHE_SparseBCF), 
                   SE_SHBCF = apply(Test_PEHE_SparseBCF, 2, function(x) MC_se(x, B)),
                   DGP = c("Not Sparse", "Sparse PS", "Sparse Mu", "Sparse Tau"))

df_tr = reshape(data = df_tr, 
                varying = list(c(1, 3), c(2, 4)), timevar = "Model", 
                v.names = c("PEHE", "SE"), direction = "long",
                times = c("BCF", "SH-BCF"), idvar = "DGP", ids = DGP)

df_te = reshape(data = df_te, 
                varying = list(c(1, 3), c(2, 4)), timevar = "Model", 
                v.names = c("PEHE", "SE"), direction = "long",
                times = c("BCF", "SH-BCF"), idvar = "DGP", ids = DGP)

df_tr; df_te


# Save
write.csv(Train_PEHE_BCF,
          paste0(curr_dir, "/Results/Single_BCF_Train_SparseDGPs_NEW_NEW.csv"))

write.csv(Test_PEHE_BCF,
          paste0(curr_dir, "/Results/Single_BCF_Test_SparseDGPs_NEW_NEW.csv"))

write.csv(Train_PEHE_SparseBCF,
          paste0(curr_dir, "/Results/Single_SHBCF_Train_SparseDGPs_NEW_NEW.csv"))

write.csv(Test_PEHE_SparseBCF,
          paste0(curr_dir, "/Results/Single_SHBCF_Test_SparseDGPs_NEW_NEW.csv"))


write.csv(df_tr,
          paste0(curr_dir, "/Results/PEHE_Train_SparseDGPs_NEW_NEW.csv"))

write.csv(df_te,
          paste0(curr_dir, "/Results/PEHE_Test_SparseDGPs_NEW_NEW.csv"))


# PLOTTING ----------------

