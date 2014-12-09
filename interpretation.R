#setwd("/accounts/class/s215a/s215a-8/Documents/final")
setwd("~/Documents/STAT215A/final_project")
library(glmnet)
library(doMC)
library(doParallel)
source("helperFunctions.R")
load("fMRIdata.RData")

###################
# Interpretation of the models
###################

##### load data
obs <- fit_feat
resp <- resp_dat
# obs and resp are divided into a training set (80%) and a validation set (20%)
select.tr <- as.logical(sample(c(1,0), size=nrow(obs), replace=TRUE, prob = c(4/5, 1/5)))
obs.train <- obs[select.tr, ]
resp.train <- resp[select.tr, ]
obs.val <- obs[!select.tr, ]
resp.val <- resp[!select.tr, ]

#### fit model
fit.5 <- cv.glmnet(x = obs.train, y = resp.train[ , 1], family = "gaussian", 
                   nlambda = 100, parallel = TRUE, type.measure = "mse",
                   nfolds = 10, alpha = 0.5)
beta <- as.numeric(coef(fit.5))[2:length(as.numeric(coef(fit.5)))]
beta[beta != 0] <- 1
beta.matrix <- matrix(rep(beta, nrow(obs)), ncol = ncol(obs), nrow = nrow(obs), 
                      byrow = TRUE)
transfer.matrix <- read.csv("real_wav.csv")
(obs * beta.matrix)

wav.best.non0.voxel1 <- ReadRealBasisFunction(28) +
  ReadRealBasisFunction(51) + ReadRealBasisFunction(253) + 
  ReadRealBasisFunction(389) + ReadRealBasisFunction(1104) + 
  ReadRealBasisFunction(1145) + ReadRealBasisFunction(1990) + 
  ReadRealBasisFunction(2044) + ReadRealBasisFunction(3891) +
  ReadRealBasisFunction(6076) + ReadRealBasisFunction(6178) +
  ReadRealBasisFunction(6382) + ReadRealBasisFunction(8804) +
  ReadRealBasisFunction(8900) + ReadRealBasisFunction(10586) +
  ReadRealBasisFunction(10615)
image(wav.best.non0.voxel1)


wav.best.non0.voxel2 <- ReadRealBasisFunction(1) + ReadRealBasisFunction(29) +
  ReadRealBasisFunction(42) + ReadRealBasisFunction(428) +
  ReadRealBasisFunction(1145) + ReadRealBasisFunction(4905) + 
  ReadRealBasisFunction(5163) + ReadRealBasisFunction(6076) + 
  ReadRealBasisFunction(6188) + ReadRealBasisFunction(6964) +
  ReadRealBasisFunction(7060) + ReadRealBasisFunction(8492) +
  ReadRealBasisFunction(9516) + ReadRealBasisFunction(10586) +
  ReadRealBasisFunction(9990)
image(wav.best.non0.voxel2)


