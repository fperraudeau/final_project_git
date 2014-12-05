setwd("~/Documents/STAT215A/final_project")
library(glmnet)
library(doMC)
library(doParallel)
library(msgps)
source("helperFunctions.R")
load("fMRIdata.RData")

###################
# GLM with Lasso, Ridge and elasticnet regularization
###################

##### Cross-validation
obs <- fit_feat
resp <- resp_dat
# obs and resp are divided into a training set (80%) and a validation set (20%)
select.tr <- as.logical(sample(c(1,0), size=nrow(obs), replace=TRUE, prob = c(4/5, 1/5)))
obs.train <- obs[select.tr, ]
resp.train <- resp[select.tr, ]
obs.val <- obs[!select.tr, ]
resp.val <- resp[!select.tr, ]

# alpha to be tested. alpha = 0 ridge regularization, 
# alpha = 1 lasso regularization
#alpha <- c(0, 0.25, 0.5, 0.75, 1)
#cv.result <- cvResult(obs.train, resp.train, obs.val, resp.val, alpha)
#write.table(cv.result, "CVresult.csv",  sep = ",", 
#            row.names = TRUE, col.names = TRUE)


#### AICc
# from ridge to lasso
alpha <- c(0.99, 0.75, 0.5, 0.25, 0)
AICc.result <- InformationCriterion(obs.train, resp.train, obs.val, resp.val, alpha, 
                     criterion = "AICc")
write.table(AICc.result, "AICcResult.csv",  sep = ",", 
            row.names = TRUE, col.names = TRUE)


#### BIC
BIC.result <- InformationCriterion(obs.train, resp.train, obs.val, resp.val, alpha, 
                                    criterion = "BIC")
write.table(BIC.result, "BICresult.csv",  sep = ",", 
            row.names = TRUE, col.names = TRUE)


#### To test on small dataset
obs.train <- obs.train[1:200, 1:500]
resp.train <- resp.train[1:200, 1:10]

obs.val <- obs.val[1:50, 1:500]
resp.val <- resp.val[1:50, 1:10]


#### ES-CV
# sequence of lambda to test, allows alignment of the lambda in the folds
# It indexes the solution path
lambda <- seq(1, 0, length.out = 100)
# number of folds
V <- 10
# from ridge to lasso
alpha <- c(0, 0.25, 0.5, 0.75, 1)
cvESCV.result <- cvESCV(lambda, V, alpha)
write.table(cvESCV.result, "cvESCV.csv",  sep = ",", 
            row.names = TRUE, col.names = TRUE)



















































