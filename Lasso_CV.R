setwd("~/Documents/STAT215A/final project")
library(glmnet)
library(doMC)
library(doParallel)
source("helperFunctions.R")
load("fMRIdata.RData")
ls()


ShrinkageCV <- function(obs, resp, nfolds, foldid, alpha){
  registerDoMC(cores=6)
  cvfit <- cv.glmnet(x = obs, y = resp, family = "mgaussian", nlambda = 100, 
                     parallel = TRUE, type.measure = "mse", nfolds = nfolds,
                     foldid = foldid, alpha = alpha)
  return(cvfit) 
}


ShrinkageResult <- function(obs.train, resp.train, obs.val, resp.val, alpha){
  ncores <- 6
  registerDoParallel(ncores)
  
  nfolds <- 10
  foldid <- sample(1:nfolds, size=nrow(resp.train), replace=TRUE)
  result <- vector()
  
  result.loop <- foreach(i = 1:length(alpha)) %dopar% {
    cvfit <- ShrinkageCV(obs.train, resp.train, nfolds = nfolds, foldid = foldid, 
                         alpha = alpha[i])
    lambda.1se <- cvfit$lambda.1se
    index.lambda.1se <- which(cvfit$lambda %in% cvfit$lambda.1se)
    pred.err <- c(cvfit$cvm[index.lambda.1se], cvfit$cvsd[index.lambda.1se]) 
    val.err <- (1 / nrow(obs.val)) * sum(
      (predict(cvfit, newx = obs.val, s = "lambda.1se")[, ,1] - resp.val)^2)
    result <- cbind(result, 
                    c(lambda.1se, pred.err[1], pred.err[2], val.err))
    return(result)
  }
  result <- as.data.frame(result.loop)
  result <- round(result, 2)
  colnames(result) <- c("Ridge", "alpha.25", "alpha.5", "alpha.75",
                        "Lasso")
  rownames(result) <- c("Lambda 1se", "Mean CV Error", "SD CV Error", "Prediction Error")
  return(as.data.frame(result))
}

obs <- fit_feat
resp <- resp_dat
select.tr <- as.logical(sample(c(1,0), size=nrow(obs), replace=TRUE, prob = c(4/5, 1/5)))
obs.train <- obs[select.tr, ]
resp.train <- resp[select.tr, ]
obs.val <- obs[!select.tr, ]
resp.val <- resp[!select.tr, ]


alpha <- c(0, 0.25, 0.5, 0.75, 1)
shrinkageResult <- ShrinkageResult(obs.train, resp.train, obs.val, resp.val, alpha)
write.table(shrinkageResult, "shrinkageResult.csv",  sep = ",", 
            row.names = TRUE, col.names = TRUE)








