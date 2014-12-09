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

#### Diagnostics
registerDoMC(cores = 1)

# fit the model again to get the beta and lambda again. We should find the same 
# beta and lambda as when we tried the different models
fit.5 <- cv.glmnet(x = obs.train, y = resp.train[ , 1], family = "gaussian", 
                   nlambda = 100, parallel = TRUE, type.measure = "mse",
                   nfolds = 10, alpha = 0.5)
fit.lasso <- cv.glmnet(x = obs.train, y = resp.train[ , 1], family = "gaussian", 
                       nlambda = 100, parallel = TRUE, type.measure = "mse",
                       nfolds = 10, alpha = 1)
# check that we found the same lambda (and so the same beta) as previously
plot(fit.5)
plot(fit.lasso)

# correlation on the entire data set
cor(predict(fit.5, newx = obs.val, s = "lambda.1se")[, 1], resp.val[ , 1])
cor(predict(fit.lasso, newx = obs.val, s = "lambda.1se")[ , 1], resp.val[ , 1])

# plot of the residuals on the entire data set
res.5 <- predict(fit.5, newx = obs, s = "lambda.1se")[ , 1] - resp[ , 1]
plot(predict(fit.5, newx = obs, s = "lambda.1se")[ , 1], res.5, 
     ylab = "Residuals", xlab = "Fitted values",
     col=ifelse(abs(res.5) >= 3, "red", "black"))
abline(h = 0, col = "blue", lwd = 2)
abline(h = 3, col = "red", lwd = 1)
abline(h = -3, col = "red", lwd = 1)

which(abs(res.5) > 3)

im1282 <- ReadImage(1282)
image(im1282, col=gray((1:500) / 501))
wav1282 <- ReadRealBasisFunction(1282)
image(wav1282)

im1313 <- ReadImage(1313)
image(im1313, col=gray((1:500) / 501))
wav1313 <- ReadRealBasisFunction(1313)
image(wav1313)


# stability of the prediction and models
# bootstrap the (resp, obs) and compute the correlation again
stability.model.5 <- stabModels(fit.5)
boxplot(stability.model.5)
mean(stability.model.5)
sd(stability.model.5)

stability.model.lasso <- stabModels(fit.lasso)
boxplot(stability.model.lasso)
mean(stability.model.lasso)
sd(stability.model.lasso)


# interprete the features
ncores <- 2
registerDoParallel(ncores)
fit.5.20 <- cv.glmnet(x = obs.train, y = resp.train[, 1:2], family = "mgaussian", 
                      nlambda = 100, parallel = TRUE, type.measure = "mse",
                      nfolds = 10, alpha = 0.5)

fit.lasso.20 <- cv.glmnet(x = obs.train, y = resp.train[, 1:2], family = "mgaussian", 
                          nlambda = 100, parallel = TRUE, type.measure = "mse",
                          nfolds = 10, alpha = 1)

# alpha = 0.5
coef.5 <- NULL
for (i in 1:ncol(resp.train)){
  coef.5 <- c(coef.5, as.matrix(coef(fit.5.20)[[i]])[,1])
}
# nrow = ncol(obs) + 1 because of the intercept
coef.matrix.5 <- matrix(coef.5, nrow = ncol(obs) + 1, ncol = 20)
coef.matrix.5[coef.matrix.5 != 0] <- 1

# lasso
coef.lasso <- NULL
for (i in 1:ncol(resp.train)){
  coef.lasso <- c(coef.lasso, as.matrix(coef(fit.lasso.20)[[i]])[,1])
}
coef.matrix.lasso <- matrix(coef.lasso, nrow = ncol(obs) + 1, ncol = 20)
coef.matrix.lasso[coef.matrix.lasso != 0] <- 1


coef.matrix <- coef.matrix.5 + coef.matrix.lasso
common.coef <- NULL
common.coef <- apply(coef.matrix, 2, function(i){which(i == 2)})
lapply(common.coef, write, "common.coef.txt", append = TRUE, ncolumns = 100000)
##### After SCF computation
both.models <- read.csv("common.coef.txt")
both.models <- both.models[510:616, 1]
# 107 features (1%)
length(both.models)/ncol(obs)



#### stability of the features
stability.features <- list()
stability.features[[1]] <- stabFeatures(0.5, 1)
stability.features[[1]] <- Reduce("+", stability.features[[1]])
stability.features[[2]] <- stabFeatures(0.5, 2)
stability.features[[2]] <- Reduce("+", stability.features[[2]])
stability.features[[3]] <- stabFeatures(1, 1)
stability.features[[3]] <- Reduce("+", stability.features[[3]])
stability.features[[4]] <- stabFeatures(1, 2)
stability.features[[4]] <- Reduce("+", stability.features[[4]])
lapply(stability.features, write, "stability.features.txt", 
       append = TRUE, ncolumns = 100000)
# after SCF computation
boot.feat <- read.csv("stability.features.txt")
boot.feat.1 <- as.vector(boot.feat[1, 1])
boot.feat.1 <- as.integer(strsplit(boot.feat.1, " ")[[1]])
boot.feat.2 <- as.vector(boot.feat[2, 1])
boot.feat.2 <- as.integer(strsplit(boot.feat.2, " ")[[1]])
boot.feat.3 <- as.vector(boot.feat[3, 1])
boot.feat.3 <- as.integer(strsplit(boot.feat.3, " ")[[1]])
boot.feat.4 <- as.vector(boot.feat[4, 1])
boot.feat.4 <- as.integer(strsplit(boot.feat.4, " ")[[1]])
which(boot.feat.3 > 700 & boot.feat.1 > 700 & boot.feat.2 >700 & boot.feat.4 > 700)
which(boot.feat.3 > 700 & boot.feat.1 > 700)
which(boot.feat.2 > 500 & boot.feat.4 > 500)











