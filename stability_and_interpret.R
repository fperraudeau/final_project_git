#setwd("/accounts/class/s215a/s215a-8/Documents/final")
setwd("~/Documents/STAT215A/final_project")
library(glmnet)
library(doMC)
library(doParallel)
source("helperFunctions.R")
load("fMRIdata.RData")

###################
# Part 2 - Interpretation of the models
###################

##### load data
obs <- fit_feat
resp <- resp_dat
# obs and resp are divided into a training set (80%) and a validation set (20%)
select.tr <- as.logical(sample(c(1,0), size=nrow(obs), replace=TRUE, 
                               prob = c(4/5, 1/5)))
obs.train <- obs[select.tr, ]
resp.train <- resp[select.tr, ]
obs.val <- obs[!select.tr, ]
resp.val <- resp[!select.tr, ]


####################
# Diagnostics
####################
registerDoMC(cores = 8)

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

# correlation on the entire data set - to check if I got the same correlation as 
# previously
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

# outliers
which(abs(res.5) > 3)
# picture 1282 has residual up to 3
im1282 <- ReadImage(1282)
image(im1282, col=gray((1:500) / 501))
wav1282 <- ReadRealBasisFunction(1282)
image(wav1282)
# picture 1313 has residual up to 3
im1313 <- ReadImage(1313)
image(im1313, col=gray((1:500) / 501))
wav1313 <- ReadRealBasisFunction(1313)
image(wav1313)


# stability of the prediction and models
# bootstrap the (resp, obs) and compute the correlation again for 1,000 bootsrap
# elastic net
stability.model.5 <- stabModels(fit.5)
boxplot(stability.model.5)
# lasso
stability.model.lasso <- stabModels(fit.lasso)
boxplot(stability.model.lasso)


#######################
# how do voxels respond to images ?
#######################
#### fit models selected for voxels 1 and 2
ncores <- 8
registerDoParallel(ncores)
fit.5.20 <- cv.glmnet(x = obs.train, y = resp.train[, 1:2], family = "mgaussian", 
                      nlambda = 100, parallel = TRUE, type.measure = "mse",
                      nfolds = 10, alpha = 0.5)

fit.lasso.20 <- cv.glmnet(x = obs.train, y = resp.train[, 1:2], family = "mgaussian", 
                          nlambda = 100, parallel = TRUE, type.measure = "mse",
                          nfolds = 10, alpha = 1)




#### Non-null coefficients for both of the models
# Run on SCF
# alpha = 0.5
coef.5 <- NULL
for (i in 1:ncol(resp.train)){
  coef.5 <- c(coef.5, as.matrix(coef(fit.5.20)[[i]])[,1])
}
# nrow = ncol(obs) + 1 because of the intercept
coef.matrix.5 <- matrix(coef.5, nrow = ncol(obs) + 1, ncol = 20)
# coef.matrix.5 has 0 if null coefficient and 1 if non null coefficient
coef.matrix.5[coef.matrix.5 != 0] <- 1

# lasso
coef.lasso <- NULL
for (i in 1:ncol(resp.train)){
  coef.lasso <- c(coef.lasso, as.matrix(coef(fit.lasso.20)[[i]])[,1])
}
# coef.matrix.lasso has 0 if null coefficient and 1 if non null coefficient
coef.matrix.lasso <- matrix(coef.lasso, nrow = ncol(obs) + 1, ncol = 20)
coef.matrix.lasso[coef.matrix.lasso != 0] <- 1

# coef.matrix has 0 if null coefficient for both of the models
# 1 if non null coefficient for one model and null coeff for the other
# and 2 if non null coeff for both of the models
coef.matrix <- coef.matrix.5 + coef.matrix.lasso
common.coef <- NULL
common.coef <- apply(coef.matrix, 2, function(i){which(i == 2)})
# write coef.matrix
lapply(common.coef, write, "common.coef.txt", append = TRUE, ncolumns = 100000)
# After SCF computation, load coef.matrix
both.models <- read.csv("common.coef.txt")
# 107 features non null common features
both.models <- both.models[510:616, 1]
# 107 features (1%)
length(both.models)/ncol(obs)





#### stability of the features
# 1,000 bootstrap samples for the two models and voxel 1 and 2
# Count the number of times over 1,000 where the features have non null
# coefficients. See description of the function stabFeatures in my_library.R
# Run on SCF
stability.features <- list()
stability.features[[1]] <- stabFeatures(0.5, 1)
stability.features[[1]] <- Reduce("+", stability.features[[1]])
stability.features[[2]] <- stabFeatures(0.5, 2)
stability.features[[2]] <- Reduce("+", stability.features[[2]])
stability.features[[3]] <- stabFeatures(1, 1)
stability.features[[3]] <- Reduce("+", stability.features[[3]])
stability.features[[4]] <- stabFeatures(1, 2)
stability.features[[4]] <- Reduce("+", stability.features[[4]])
# write results
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

# real part of the Gabor wavelet plotted as an image
# voxel 1 
which(boot.feat.3 > 765 & boot.feat.1 > 765)
wav.best.non0.voxel1 <- ReadRealBasisFunction(51) + 
  ReadRealBasisFunction(253) + 
  ReadRealBasisFunction(389) +  
  ReadRealBasisFunction(1145) + 
  ReadRealBasisFunction(2044) + 
  ReadRealBasisFunction(6076) + 
  ReadRealBasisFunction(6382) + 
  ReadRealBasisFunction(10586) +
  ReadRealBasisFunction(10615)
image(wav.best.non0.voxel1)
# voxel 2
which(boot.feat.2 > 610 & boot.feat.4 > 610)
wav.best.non0.voxel2.600 <- ReadRealBasisFunction(1) +
  ReadRealBasisFunction(428) +
  ReadRealBasisFunction(1145) + ReadRealBasisFunction(6076) + 
  ReadRealBasisFunction(6188) + ReadRealBasisFunction(6964) +
  ReadRealBasisFunction(9516) + ReadRealBasisFunction(10586) 
image(wav.best.non0.voxel2.600)


##### clustering the voxels
coef.clust <- NULL
registerDoMC(cores = 8)
for (i in 1:ncol(obs)){
  fit.5 <- cv.glmnet(x = obs, y = resp[ , i], family = "gaussian", 
                     nlambda = 100, parallel = TRUE, type.measure = "mse",
                     nfolds = 10, alpha = 0.5)
  coef.clust <- cbind(coef.clust, coef(fit.5)[, 1])
}
write.table(coef.clust, "coef.clust.csv",  sep = ",", 
            row.names = TRUE, col.names = TRUE)


# kmeans to cluster the voxels according to the coefficient
coef.kmean <- kmeans(t(coef.clust), centers = 5)
voxel.clust <- coef.kmean$cluster
write.table(voxel.clust, "voxel.clust.csv",  sep = ",", 
            row.names = TRUE, col.names = TRUE)


# Plot the physical locations of the voxels.
voxel.clust <- read.csv("voxel.clust5.csv")
voxel.locs <- data.frame(loc_dat)
rgl.spheres(voxel.locs$X1, voxel.locs$X2, voxel.locs$X3,
            color=voxel.clust[, "x"], radius=0.3)







