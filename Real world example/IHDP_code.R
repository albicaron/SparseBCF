##################################################
# Effects of early educational and family support on cognitive function of low birth weight preterm infants  
# original paper DOI: https://doi.org/10.1056/NEJM199610103351501
##################################################

rm(list = ls())


### LIBRARIES
library(tidyverse)
library(SparseBCF) # OOB Sp. Bayesian Causal Forests
library(BART) # Main package including all the version of BART
library(nnet)
library(rpart)
library(rpart.plot)

# Load data from Hill (2011)
script.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
load(paste0(script.dir, "/example.data"))

# Filter out duplicate variables
df = ihdp[, -sort(c(grep(".1", colnames(ihdp)), grep("F", colnames(ihdp))))[-c(1, 11)]]
df = df[, -2]
df = na.omit(df)

summary(ihdp)
summary(df)


# Relabel some variables
df = df %>%
  rename(bw_2000 = bwg) %>%
  mutate(othstudy = ifelse(othstudy == 2, 1, 0),
         drugs = ifelse(drugs == 2, 1, 0),
         language = ifelse(language == 1, 1, 0),
         momed4F = as.integer(momed4F)) %>%
  select(- c(mom.lths, mom.hs, mom.coll, mom.scoll)) # for some variable we keep factors instead of one-hot encoding 

summary(df)


################# ANALYSIS  --------------------
Y = df$iqsb.36
Z = df$treat
X = as.matrix(df[, -c(1, 2)])
P = ncol(X)

summary(df)

set.seed(1234)


### PScore Model - 1 hidden layer neural net
PS_nn <- nnet(x = X, y = Z, size = 20, maxit = 2000, 
              decay = 0.01, trace=FALSE, abstol = 1.0e-8) 
PS_est = PS_nn$fitted.values

# Check NN performance 
MLmetrics::AUC(PS_est, Z)


### Run normal BCF
bcf <- SparseBCF(y = Y, z = Z, x_control = X, 
                 pihat = PS_est, 
                 OOB = F,
                 sparse = F,
                 update_interval = 1000,
                 nburn = 10000, nsim = 5000)


### Run Sparse BCF
SPbcf <- SparseBCF(y = Y, z = Z, x_control = X, 
                   pihat = PS_est, 
                   OOB = F,
                   sparse = T,
                   update_interval = 1000,
                   nburn = 10000, nsim = 5000)

tau_BCF = colMeans(bcf$tau)
tau_SPBCF = colMeans(SPbcf$tau)

mean(tau_BCF); mean(tau_SPBCF)
sd(tau_BCF); sd(tau_SPBCF)

hist(tau_BCF)
hist(tau_SPBCF)

barplot(colMeans(SPbcf$varprb_mu))
barplot(colMeans(SPbcf$varprb_tau))


########### Plottin ##############
# Get the 10 PS percentiles
p = quantile(PS_est, probs = seq(0, 1, 0.1))
index = c(NA)

for (i in 1:length(p)) {
  index[i] = which(abs(PS_est - p[i]) == (min(abs(PS_est - p[i]))))
}


# melt dataframe
BCF_est = as.data.frame(bcf$tau[, index])
SPBCF_est = as.data.frame(SPbcf$tau[, index])

colnames(BCF_est) = names(p)
colnames(SPBCF_est) = names(p)

BCF_est = reshape2::melt(BCF_est)
SPBCF_est = reshape2::melt(SPBCF_est)

colnames(BCF_est) = c("PS_cent", "BCF")
colnames(SPBCF_est) = c("PS_cent", "SPBCF")

df = data.frame(BCF_est, SPBCF_est)
df = df[-3]

df_final = reshape2::melt(df)
colnames(df_final) = c("PS_cent", "Model", "CATE")


# Joy plots
library(ggridges)

df_ = df_final[df_final$Model == "SPBCF", ]

CATEplot <-
  ggplot(df_, aes(y = PS_cent, x = CATE, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE, scale = 1.5, size = 0.75, rel_min_height = 0.025) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  # scale_fill_viridis_c(name = "Tail Probability", begin = 0.1, direction = -1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(limits = c(-10, 30), breaks = seq(-10, 30, 5)) + xlab("CATE") +
  theme_ridges(center_axis_labels = TRUE) + ylab("Propensity Percentile") +
  theme(legend.position = "none", text = element_text(size=13.5))

CATEplot



########### Other plots
TAU <- 
  ggplot(data = data.frame(tau_SPBCF)) + 
  geom_density(aes(x = tau_SPBCF, y = ..count../sum(..count..), fill = "SP-BCF"), alpha = 0.6, size = 0.8) + 
  geom_vline(xintercept = mean(tau_SPBCF), linetype = "dashed", size = 0.8) +
  geom_vline(xintercept = 0, linetype = "dotted", size = 0.5) +
  scale_fill_manual(name = "", values = c("SP-BCF" = "blue1")) +
  theme_minimal() + scale_x_continuous(breaks = round(seq(-2, 20, 2), 2), limits = c(-2, 20)) +
  ylim(c(0, 0.01)) + ylab("") + xlab("CATE") +
  theme(text = element_text(size=13.5))

TAU


Tau_25 <- 
  ggplot(data = data.frame(X = 1:P, Tau = colMeans(SPbcf$varprb_tau)), aes(X, Tau)) + 
  geom_point(size = 2, color = "blue") + geom_segment(aes(x=X, xend=X, y=0, yend=Tau, color = "SP-BCF"), size = 0.8) +
  geom_hline(aes(yintercept = 1/P, color = "BCF"), linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, color = "black") +
  theme_minimal() + scale_y_continuous(breaks = seq(0, 0.2, 0.01)) + 
  scale_x_continuous(limits = c(1, P), minor_breaks = seq(1, P, 2), breaks = seq(1, P, 2)) +
  theme(text = element_text(size=14.5)) +
  scale_color_manual(name = "", values = c("SP-BCF" = "blue", "BCF" = "gray41")) +
  ylab(expression("Splitting Probability")) + xlab("j-th predictor") + 
  guides(color = guide_legend(override.aes = list(size = 1)))

Tau_25



# Saving
ggsave("YOUR DIRECTORY", 
       ggpubr::ggarrange(TAU, Tau_25), 
       width = 32, height = 12, units = "cm")

# Saving
ggsave("YOUR DIRECTORY", 
       ggpubr::ggarrange(CATEplot, Tau_25), 
       width = 32, height = 12, units = "cm")


######## Plot partition tree

# Change relevantvariables labels ex-post
colnames(X)[which(colnames(X) %in% c("parity", "bw_2000", "momwhite", "momed4F"))] = 
  c("Num of children", "Birth weight under 2kg", "White mom", "Mom's educ level")

mytree <- rpart(
  tau_SPBCF ~ ., 
  data =  cbind.data.frame(tau_SPBCF, X), 
  control = list(maxdepth = 3)
)

# pdf("YOUR DIRECTORY", 
#     width = 12, height = 8)
rpart.plot(mytree, type = 2, extra = 101, clip.right.labs = FALSE, 
           box.palette = "Blues", # color scheme
           branch.lty = 3, # dotted branch lines
           shadow.col = "gray",
           branch.lwd = 2,
           tweak = 1.1,
           branch = 1, under = TRUE,  yesno = 2)
# dev.off()

