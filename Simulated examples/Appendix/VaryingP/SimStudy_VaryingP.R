###################################################
# Some other methods to compare on Simulated Data #
###################################################

rm(list = ls())


### LIBRARIES
library(tidyverse)
library(SparseBCF) # OOB Bayesian Causal Forests
library(BART) # Main package including all the version of BART

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
  X[, 1:ceiling(P/2)] = qnorm(unif[, 1:ceiling(P/2)])
  X[, (ceiling(P/2)+1):P] = qbinom(unif[, (ceiling(P/2)+1):P], 1, 0.3)
  
  return(X)
  
}


### INFO
N = 250
P = c(5, 10, 50, 100, 150)
B = 200


# Store simulations true and estimated quantities
Train_PEHE_BCF = matrix(NA, B, length(P)); Test_PEHE_BCF = matrix(NA, B, length(P))
Train_PEHE_SparseBCF = matrix(NA, B, length(P)); Test_PEHE_SparseBCF = matrix(NA, B, length(P))

for (j in 1:length(P)) {
  for (i in 1:B) {
    
    
    gc()
    
    cat("\n\n\n\n-------- Iteration", i, "--------\n\n\n\n")
    
    
    if(i<=500){set.seed(100 + i*5)}
    if(i>500){set.seed(5000 + i*5)}
    
    
    # Simulate covariates
    X <- get_features(N, P[j])
    
    
    # Simulate Pscore and Z
    first_bin = (ceiling(P[j]/2)+1)
    
    und_lin = -0.5 + 0.2*X[, 1] + 0.1*X[, 2] + 0.4*X[, first_bin] + runif(N)/10
    pscore <- pnorm( und_lin )
    
    Z = rbinom(N, 1, pscore)
    
    
    ############## RESPONSE SURFACES
    mu = 3 + 1.5*sin(pi*X[, 1]) + 0.5*(X[, 2] - 0.5)^2 + 1.5*(2-abs(X[, 3])) + 1.5*X[, 4]*(X[, first_bin] + 1)
    ITE = 0.1 +  1*abs(X[, 1] - 1)*(X[, first_bin] + 2)
    
    sigma = sd(mu)/2
    Y = mu + ITE*Z + rnorm(N, 0, sigma)
    
    
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
                   P = P)

df_te = data.frame(PEHE_BCF = colMeans(Test_PEHE_BCF), 
                   SE_BCF = apply(Test_PEHE_BCF, 2, function(x) MC_se(x, B)),
                   PEHE_SHBCF = colMeans(Test_PEHE_SparseBCF), 
                   SE_SHBCF = apply(Test_PEHE_SparseBCF, 2, function(x) MC_se(x, B)),
                   P = P)

df_tr = reshape(data = df_tr, 
                varying = list(c(1, 3), c(2, 4)), timevar = "Model", 
                v.names = c("PEHE", "SE"), direction = "long",
                times = c("BCF", "SH-BCF"), idvar = "P", ids = P)

df_te = reshape(data = df_te, 
                varying = list(c(1, 3), c(2, 4)), timevar = "Model", 
                v.names = c("PEHE", "SE"), direction = "long",
                times = c("BCF", "SH-BCF"), idvar = "P", ids = P)

df_tr; df_te


# Save
curr_dir <- dirname(rstudioapi::getSourceEditorContext()$path)


write.csv(Train_PEHE_BCF,
          paste0(curr_dir, "/Results/Single_BCF_Train_VaryingP_NEW_NEW.csv"))

write.csv(Test_PEHE_BCF,
          paste0(curr_dir, "/Results/Single_BCF_Test_VaryingP_NEW_NEW.csv"))

write.csv(Train_PEHE_SparseBCF,
          paste0(curr_dir, "/Results/Single_SHBCF_Train_VaryingP_NEW_NEW.csv"))

write.csv(Test_PEHE_SparseBCF,
          paste0(curr_dir, "/Results/Single_SHBCF_Test_VaryingP_NEW_NEW.csv"))


write.csv(df_tr,
          paste0(curr_dir, "/Results/PEHE_Train_VaryingP_NEW_NEW.csv"))

write.csv(df_te,
          paste0(curr_dir, "/Results/PEHE_Test_VaryingP_NEW_NEW.csv"))



# PLOTTING ----------------
colnames(Train_PEHE_BCF) <- P; colnames(Test_PEHE_BCF) <- P
colnames(Train_PEHE_SparseBCF) <- P; colnames(Test_PEHE_SparseBCF) <- P

all_tr_BCF = reshape2::melt(Train_PEHE_BCF)
all_tr_BCF$Var1 = "BCF"
all_tr_BCF$Var2 = as.factor(all_tr_BCF$Var2)
colnames(all_tr_BCF) = c("Model", "P", "PEHE")

all_te_BCF = reshape2::melt(Test_PEHE_BCF)
all_te_BCF$Var1 = "BCF"
all_te_BCF$Var2 = as.factor(all_te_BCF$Var2)
colnames(all_te_BCF) = c("Model", "P", "PEHE")

all_tr_SparseBCF = reshape2::melt(Train_PEHE_SparseBCF)
all_tr_SparseBCF$Var1 = "SparseBCF"
all_tr_SparseBCF$Var2 = as.factor(all_tr_SparseBCF$Var2)
colnames(all_tr_SparseBCF) = c("Model", "P", "PEHE")

all_te_SparseBCF = reshape2::melt(Test_PEHE_SparseBCF)
all_te_SparseBCF$Var1 = "SparseBCF"
all_te_SparseBCF$Var2 = as.factor(all_te_SparseBCF$Var2)
colnames(all_te_SparseBCF) = c("Model", "P", "PEHE")


all_tr = rbind(all_tr_BCF, all_tr_SparseBCF)
all_te = rbind(all_te_BCF, all_te_SparseBCF)


library(ggridges)
library(ggar)

p1 =
  ggplot(all_tr, aes(x=PEHE, y=P, fill=Model, linetype=Model)) +
  geom_density_ridges(scale=0.8, rel_min_height=.025,
                      alpha = 0.5, point_alpha=1, size = 0.8) + 
  theme_minimal() + coord_flip() + ylab("Num of predictors") + 
  xlab(expression('Train ' ~ sqrt(PEHE))) + scale_x_continuous(breaks = seq(0, 2.5, 0.5)) +
  scale_fill_manual(values = c("#4738FA", "#D64B00")) + xlim(c(0.2, 2.5)) +
  theme(text = element_text(size=13)) +
  theme(legend.key.size = unit(0.7, 'cm'), #change legend key size
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=9)) #change legend text font size
p1

p2 =
  ggplot(all_te, aes(x=PEHE, y=P, fill=Model, linetype=Model)) +
  geom_density_ridges(scale=0.8, rel_min_height=.025,
                      alpha = 0.5, point_alpha=1, size = 0.8) + 
  theme_minimal() + coord_flip() + ylab("Num of predictors") + 
  xlab(expression('Test ' ~ sqrt(PEHE))) + scale_x_continuous(breaks = seq(0, 2.5, 0.5)) +
  scale_fill_manual(values = c("#4738FA", "#D64B00")) + xlim(c(0.2, 2.5)) +
  theme(text = element_text(size=13)) +
  theme(legend.key.size = unit(0.7, 'cm'), #change legend key size
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=9)) #change legend text font size
p2


ggsave(paste0(curr_dir, "/Results/TrainTest_VaryingP.pdf"), 
       ggpubr::ggarrange(p1, p2), 
       width = 30, height = 12, units = "cm")



