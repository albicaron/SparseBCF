###########################
# Simulated Guide Example #
###########################

rm(list = ls())

# Libraries
library(tidyverse)
library(BART) # Main package including all the version of BART
library(SparseBCF)
library(rlearner)
library(plotly)

# Functions
PEHE <- function(x,y) sqrt(mean((x-y)^2))

# Options
P = 5
N = 1000
set.seed(1234)



# Simulate correlated covariates -----------------------------
mysigma = matrix(1, P, P)

for (i in 1:P) {
  for (j in 1:P) {
    mysigma[i, j] = 0.6^abs(i - j) + ifelse(i == j, 0, 0.1)
  }
}


X = MASS::mvrnorm(N, mu = rep(0, P), mysigma)
z <- rbinom(N, 1, pnorm(-0.5 + 0.4*X[, 1]) )

summary(pnorm(-0.5 + 0.4*X[, 1]))
cor(X)

# Simulate POs and Y
mu = 3 + X[, 1]
ITE = 0.5 + 0.5*X[, 2]^2

Y = mu + ITE*z + rnorm(N, 0, 1)

summary(Y)
summary(mu); summary(ITE)


# S-DART -----------------------------------------
SDART <-
  wbart(cbind(X, z), Y, sparse = T, 
        nskip = 1000, ndpost = 2000)

SDART_cnt <- colMeans( predict(SDART, cbind(X, ifelse(z == 1, 0, 1))) )

SDART0 = rep(NA, N); SDART1 = rep(NA, N)
SDART0[z==0] = SDART$yhat.train.mean[z==0]
SDART0[z!=0] = SDART_cnt[z!=0]

SDART1[z==1] = SDART$yhat.train.mean[z==1]
SDART1[z!=1] = SDART_cnt[z!=1]

tau_SDART = SDART1 - SDART0
  

# T-DART ----------------------------------------
TDART0 <-
  wbart(X[z==0, ], Y[z==0], sparse = T,
        nskip = 1000, ndpost = 2000)

TDART1 <-
  wbart(X[z==1, ], Y[z==1], sparse = T,
        nskip = 1000, ndpost = 2000)


Y0TDART = rep(NA, N); Y1TDART = rep(NA, N)

Y0TDART[z==0] = TDART0$yhat.train.mean
Y1TDART[z==1] = TDART1$yhat.train.mean

TDART0_fit <- colMeans(predict(TDART0, as.matrix(X)))
TDART1_fit <- colMeans(predict(TDART1, as.matrix(X)))

tau_TDART = TDART1_fit - TDART0_fit


# SDART and TDART are unable to carry out variable selection on the moderating effects
# PS estimation (we use Probit BART in this case as classifier, just to show different models are suitable)
PS_bart <- pbart(X, z, nskip = 1000, ndpost = 2000)
PS_est = PS_bart$prob.train.mean

# Sparse BCF ---------------------------------------
mysbcf <-
  SparseBCF(y = Y, z = z, x_control = X, 
            pihat = PS_est, 
            OOB = F, 
            sparse = T,
            update_interval = 2000,
            nburn = 1000, nsim = 2000)


tau_SBCF = colMeans(mysbcf$tau)


colMeans(SDART$varprob)

colMeans(TDART0$varprob)
colMeans(TDART1$varprob)

colMeans(mysbcf$varprb_mu)
colMeans(mysbcf$varprb_tau)


# PEHE
PEHE(tau_SDART, ITE)
PEHE(tau_TDART, ITE)
PEHE(tau_SBCF, ITE)

