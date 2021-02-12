######################################
# Sparse BCF splitting probabilities #
######################################

rm(list = ls())

# Libraries
library(tidyverse)
library(latex2exp)

# Load Spliting proba
# WARNING: THIS CODE TAKES RESULTS OF "AllR_Models_P25.R" and "AllR_Models_P50.R" AS INPUTS
# TO PRODUCE LOLLIPOP CHARTS
splitBCF_25 <- read.csv("YOUR DIRECTORY/Split_Final_P25.csv")
splitBCF_50 <- read.csv("YOUR DIRECTORY/Split_Final_P50.csv")



# Plot
Mu_25 <- 
  ggplot(data = splitBCF_25, aes(X, Mu)) + 
  geom_point(size = 1.7, color = "brown1") + geom_segment(aes(x=X, xend=X, y=0, yend=Mu, color = "SP-BCF"), size = 0.6) +
  geom_hline(aes(yintercept = 0.04, color = "BCF"), linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, color = "black") +
  theme_light() + scale_y_continuous(breaks = seq(0, 0.3, 0.05)) + 
  scale_x_continuous(limits = c(1, 25), minor_breaks = seq(1, 25, 2), breaks = seq(1, 25, 2)) +
  theme(text = element_text(size=14.5)) + 
  scale_color_manual(name = "", values = c("SP-BCF" = "brown1", "BCF" = "gray41")) +
  ylab(expression("s" [mu])) + xlab("j-th predictor") + guides(color = guide_legend(override.aes = list(size = 1)))

Mu_25


Tau_25 <- 
  ggplot(data = splitBCF_25, aes(X, Tau)) + 
  geom_point(size = 1.7, color = "cornflowerblue") + geom_segment(aes(x=X, xend=X, y=0, yend=Tau, color = "SP-BCF"), size = 0.6) +
  geom_hline(aes(yintercept = 0.04, color = "BCF"), linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, color = "black") +
  theme_light() + scale_y_continuous(breaks = seq(0, 0.5, 0.1)) + 
  scale_x_continuous(limits = c(1, 25), minor_breaks = seq(1, 25, 2), breaks = seq(1, 25, 2)) +
  theme(text = element_text(size=14.5)) +
  scale_color_manual(name = "", values = c("SP-BCF" = "cornflowerblue", "BCF" = "gray41")) +
  ylab(expression("s" [tau])) + xlab("j-th predictor") + guides(color = guide_legend(override.aes = list(size = 1)))

Tau_25


Mu_50 <- 
  ggplot(data = splitBCF_50, aes(X, Mu)) + 
  geom_point(size = 1.7, color = "brown1") + geom_segment(aes(x=X, xend=X, y=0, yend=Mu, color = "SP-BCF"), size = 0.6) +
  geom_hline(aes(yintercept = 0.02, color = "BCF"), linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, color = "black") +
  theme_light() + scale_y_continuous(breaks = seq(0, 0.3, 0.05)) + 
  scale_x_continuous(limits = c(1, 50), minor_breaks = seq(1, 50, 4), breaks = seq(1, 50, 4)) +
  theme(text = element_text(size=14.5)) +
  scale_color_manual(name = "", values = c("SP-BCF" = "brown1", "BCF" = "gray41")) +
  ylab(expression("s" [mu])) + xlab("j-th predictor") + guides(color = guide_legend(override.aes = list(size = 1)))

Mu_50


Tau_50 <- 
  ggplot(data = splitBCF_50, aes(X, Tau)) + 
  geom_point(size = 1.7, color = "cornflowerblue") + geom_segment(aes(x=X, xend=X, y=0, yend=Tau, color = "SP-BCF"), size = 0.6) +
  geom_hline(aes(yintercept = 0.02, color = "BCF"), linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = 0, color = "black") +
  theme_light() + scale_y_continuous(breaks = seq(0, 0.5, 0.1)) + 
  scale_x_continuous(limits = c(1, 50), minor_breaks = seq(1, 50, 4), breaks = seq(1, 50, 4)) +
  theme(text = element_text(size=14.5)) +
  scale_color_manual(name = "", values = c("SP-BCF" = "cornflowerblue", "BCF" = "gray41")) +
  ylab(expression("s" [tau])) + xlab("j-th predictor") + guides(color = guide_legend(override.aes = list(size = 1)))

Tau_50



# Save plots
ggsave("YOUR DIRECTORY/Mu_25.pdf", 
       plot = Mu_25, device = "pdf", width = 14, height = 10, units = "cm")

ggsave("YOUR DIRECTORY/Tau_25.pdf", 
       plot = Tau_25, device = "pdf", width = 14, height = 10, units = "cm")

ggsave("YOUR DIRECTORY/Mu_50.pdf", 
       plot = Mu_50, device = "pdf", width = 14, height = 10, units = "cm")

ggsave("YOUR DIRECTORY/Tau_50.pdf", 
       plot = Tau_50, device = "pdf", width = 14, height = 10, units = "cm")
