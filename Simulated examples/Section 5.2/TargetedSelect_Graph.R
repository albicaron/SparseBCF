####################################
# Targeted Selection via SparseBCF #
####################################

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
P = 10
set.seed(123)

# Simulate correlated covariates
X <- get_features(N, P)


# /////////////////////////
# 2) Targeted Selection ----------------------
# /////////////////////////
# Generate quantities, with strong targeted selection
mu_base = 2 + 0.5*sin(pi*X[, 1]) - 0.25*X[, 2]^2 + 0.75*X[, 3]*X[, 9]  
mu = 5*mu_base
# plot(X[, 1],  mu)    # Uncomment to check quantities

ITE = 1 + 2*abs(X[, 4]) + 1*X[, 10]
# plot(X[, 4], ITE)           # Uncomment to check quantities
# summary(mu); summary(ITE)    
# par(mfrow = c(1, 2))
# hist(mu); hist(ITE)

Pscore = 0.9*plogis(1.2 - mu_base + runif(N)/10)
z <- rbinom(N, 1, Pscore)
# summary(Pscore)   # Uncomment to check quantities
# hist(Pscore)
# table(z)

Y = mu + ITE*z + rnorm(N, 0, sd(ITE)/2)
# summary(Y)   # Uncomment to check quantities
# summary(mu); summary(ITE)


# Estimate Pscore (this time with 1 layers NN)  --------------------------
PS_nn <- nnet(x = X, y = z, size = 12, maxit = 2000, 
              decay = 0.01, trace=FALSE, abstol = 1.0e-8) 
PS_est = PS_nn$fitted.values

# Check NN performance 
MLmetrics::AUC(PS_est, z)


# Normal BCF
TarSel_BCF <-
  SparseBCF(y = Y, z = z, x_control = X, 
            pihat = PS_est, 
            OOB = F,
            sparse = F,
            update_interval = 1000,
            nburn = 10000, nsim = 5000)


# Sparse BCF -------------------------
TarSel_SBCP <-
  SparseBCF(y = Y, z = z, x_control = X, 
            pihat = PS_est, 
            OOB = F,
            sparse = T,
            update_interval = 1000,
            nburn = 10000, nsim = 5000)


# Informative Sparse BCF -------------------------
TarSel_IBCP <-
  SparseBCF(y = Y, z = z, x_control = X, 
            pihat = PS_est, 
            OOB = F,
            sparse = T, 
            inform_mu = T, weights_mu = c(rep(1, P), 50),
            update_interval = 1000,
            nburn = 10000, nsim = 5000)


TarSel_IBCP2 <-
  SparseBCF(y = Y, z = z, x_control = X, 
            pihat = PS_est, 
            OOB = F,
            sparse = T, 
            inform_mu = T, weights_mu = c(rep(1, P), 100),
            update_interval = 1000,
            nburn = 10000, nsim = 5000)


colMeans(TarSel_BCF$varprb_mu)[P+1]; colMeans(TarSel_SBCP$varprb_mu)[P+1];
colMeans(TarSel_IBCP$varprb_mu)[P+1]; colMeans(TarSel_IBCP2$varprb_mu)[P+1];

# Plot PS-mu fit
df = data.frame(PS = Pscore, mu = mu,
                mu_bcf = colMeans(TarSel_BCF$mu), mu_sbcf = colMeans(TarSel_SBCP$mu),
                mu_ibcf = colMeans(TarSel_IBCP$mu), mu_ibcf2 = colMeans(TarSel_IBCP2$mu))

df = reshape(data = df, varying = list(3:ncol(df)), timevar = "Type", v.names = "mu_fit", direction = "long")
df = df %>% 
  select(-id) %>% 
  mutate(Type = factor(Type, ordered = T, labels = c("BCF", "SP-BCF", "I-BCF (50)", "I-BCF (100)")))


myColors <- list("BCF" = "yellowgreen", "SP-BCF" = "#FF5000", 
                 "I-BCF (50)" = "#77BBFF", "I-BCF (100)" = "#566FFB")

myplot = 
  ggplot(df) + geom_point(aes(PS,mu, fill = "True"), alpha = 0.8) +
  scale_fill_manual(name = "", values = c("True" = "black")) + 
  geom_point(aes(PS, mu_fit, color = Type), alpha = 0.5) + theme_minimal() + scale_colour_manual(name = 'Fit:', values = myColors)+ 
  theme(text = element_text(size=15)) +
  facet_wrap(~ Type, nrow = 2, ) + ylab("Mu") + xlab("PS") + 
  guides(fill = guide_legend(override.aes = list(size=2.5)), color = guide_legend(override.aes = list(alpha=0.8, size=3)))

ggsave(filename = "YOUR DIRECTORY/MUvsPS_TarSel.pdf", 
       plot = myplot, device = "pdf", height = 8, width = 11, units = "in")


# Check Mu and ITE RMSE
PEHE(mu, colMeans(TarSel_BCF$mu))
PEHE(mu, colMeans(TarSel_SBCP$mu))
PEHE(mu, colMeans(TarSel_IBCP$mu))
PEHE(mu, colMeans(TarSel_IBCP2$mu))

PEHE(ITE, colMeans(TarSel_BCF$tau))
PEHE(ITE, colMeans(TarSel_SBCP$tau))
PEHE(ITE, colMeans(TarSel_IBCP$tau))
PEHE(ITE, colMeans(TarSel_IBCP2$tau))
