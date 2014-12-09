### This file contains all the functions used for the final project

RotateImageVector <- function(image.vec) {
  # Rotate a raw image vector into the
  # correct orientation for viewing.
  
  # There are this many rows and columns,
  # and transpose and re-ordering gives it the correct
  # orientation.
  kRows <- 128
  return(t(matrix(image.vec, nrow = kRows)[kRows:1,]))
}


ReadImage <- function(number) {
  # Read a single line from fit_stim.csv into
  # a matrix that can be plotted as an image.
  #
  # Args:
  #  - number: Which image number to load, starting
  #            with 1.
  #
  # Output:
  #  - A matrix containing the raw image pixels.
  
  # The first line contains column names, so we always
  # skip at least one line.
  img <- scan("fit_stim.csv",
              skip = number, nlines = 1, sep=",")
  
  return(RotateImageVector(img))
  
}

ReadRealBasisFunction <- function(number) {
  # Read a single column from real_wav.csv into
  # a matrix that can be plotted as an image.  This
  # uses a bash command and may not work for windows users. 
  #
  # Args:
  #  - number: Which basis function to load, starting
  #            with 1.
  #
  # Output:
  #  - A matrix containing the basis function.  NB:
  #    I am not 100% sure that the rotation is the same
  #    for the basis function as for the images.
  
  # Use bash to read a single column into a text file.
  # Note that this may not work for windows users.
  temp.file <- tempfile()
  system(paste("cat real_wav.csv | cut -d, -f ",
               number, " > ",
               temp.file))
  basis.data <- read.csv(temp.file, header=T)[[1]]
  return(RotateImageVector(basis.data))
}


crossValidation <- function(obs, resp, nfolds, foldid, alpha){
  # cross-validation to find the lambda for a GLM with lasso or elasticnet 
  # regularization.
  # To be run in parallel, 6 cores are used
  #
  # Args:
  #   - obs : matrix of the observations
  #   - resp : matrix of the responses
  #   - nfolds : integer,number of fold to be used for the cross validation
  #   - alpha : numeric, elasticnet mixing parameter. If alpha = 0, the model is 
  #   a generalized linear model with ridge regularization. If alpha = 1, the 
  #   model is a generalized linear model with lasso regularization
  # 
  # Output:
  #   - cvfit : an object of class "cv.glmnet" is returned, which is a list with
  #   the ingredients of the cross-validation fit.
  
  # use 6 cores
  registerDoMC(cores = 2)
  # cross validation with family = "mgaussian" because the response is a 
  # matrix with 20 voxels (ie 20 columns)
  cvfit <- cv.glmnet(x = obs, y = resp, family = "mgaussian", nlambda = 100, 
                     parallel = TRUE, type.measure = "mse", nfolds = nfolds,
                     foldid = foldid, alpha = alpha)
  return(cvfit) 
}


cvResult <- function(obs.train, resp.train, obs.val, resp.val, alpha){
  # cross-validation to find the lambda for a GLM with lasso or elastic net 
  # regularization.
  # To be run in parallel, 8 cores are used
  #
  # Args:
  #   - obs.train : matrix of the observations for the training set
  #   - resp.train : matrix of the responses for the training set
  #   - obs.val : matrix of the observations for the validation set
  #   - resp.val : matrix of the responses for the validation set
  #   - alpha : numeric, elasticnet mixing parameter. If alpha = 0, the model is 
  #   a generalized linear model with ridge regularization. If alpha = 1, the 
  #   model is a generalized linear model with lasso regularization
  #
  # Output:
  #   - cvfit : an object of class "cv.glmnet" is returned, which is a list with
  #   the ingredients of the cross-validation fit.
  
  # use 8 cores to run in parallel
  ncores <- 8
  registerDoParallel(ncores)
  
  # cross validation with 10 folds
  nfolds <- 10
  # foldid fixed to run cross validation with the same folds for all the alpha.
  # We want to have the same fold to compare GLM with lasso, elasticnet and ridge
  # regularization.
  foldid <- sample(1:nfolds, size=nrow(resp.train), replace=TRUE)
  result <- vector()
  
  result.loop <- foreach(i = 1:length(alpha)) %dopar% {
    # loop over the alpha to perform cross validation in parallel
    cvfit <- crossValidation(obs.train, resp.train, nfolds = nfolds, 
                             foldid = foldid, alpha = alpha[i])
    # lambda.1se is the largest value of lambda such that error is within 1 
    # standard error of the minimum lambda. 
    lambda.min <- cvfit$lambda.min
    index.lambda.min <- which(cvfit$lambda %in% cvfit$lambda.min)
    # number of non-zero coefficients at lambda.1se
    non.zero <- cvfit$nzero[index.lambda.min]
    names(non.zero) <- NULL
    # mean predicted error and standard deviation of the predicted error
    est.err <- c(cvfit$cvm[index.lambda.min], cvfit$cvsd[index.lambda.min])
    # MSE on the validation set with lambda.1se
    pred.err <- (1 / nrow(obs.val)) * sum(
      (predict(cvfit, newx = obs.val, s = "lambda.1se")[, ,1] - resp.val)^2)
    corr <- cor(predict(cvfit, newx = obs.val, s = "lambda.1se")[, ,1], resp.val)
    result <- cbind(result, 
                    c(lambda.min, est.err[1], est.err[2], pred.err, non.zero, corr))
    return(result)
  }
  result <- as.data.frame(result.loop)
  result <- round(result, 2)
  colnames(result) <- c("Ridge", "alpha.25", "alpha.5", "alpha.75",
                        "Lasso")
  rownames(result) <- c("Lambda 1se", "Mean CV Estimation Error", 
                        "SD CV Estimation Error", 
                        "Prediction Error", "Non-zero coefficient", 
                        "Correlation")
  return(as.data.frame(result))
}



InformationCriterion <- function(obs.train, resp.train, obs.val, resp.val, 
                                 alpha, criterion){
  # Return the results for AICc and BIC with the minimum lambda, the correlation
  # between the fitted and observed values and the number of non null 
  # coefficients
  # Args: 
  #   - obs.train : matrix of the observations for the training set
  #   - resp.train : matrix of the responses for the training set
  #   - obs.val : matrix of the observations for the validation set
  #   - resp.val : matrix of the responses for the validation set
  #   - alpha : numeric, elasticnet mixing parameter. If alpha = 0, the model is 
  #   a generalized linear model with ridge regularization. If alpha = 1, the 
  #   model is a generalized linear model with lasso regularization
  #   criterion : "AICc" or "BIC", the criterion used to select to model
  # Output:
  #   - a matrix with the input lambda, the correlation score, the prediction
  #    error on the validation set and the number of non null coefficient with
  #    lambda for each alpha
  lambda <- NULL
  df <- NULL
  est.err <- NULL
  pred.err <- NULL
  ave.df <- NULL
  ave.lambda <- NULL
  ave.est.err <- NULL
  ave.pred.err <- NULL
  corr <- NULL
  ave.corr <- NULL
  if (criterion == "AICc"){
    k <- 5
    df_result <- "dfaicc_result"
  }
  if (criterion == "BIC"){
    k <- 7
    df_result <- "dfbic_result"
  }
  for (i in 1:length(alpha)){
    for (j in 1:ncol(resp.train)){
      fit <- msgps(obs.train, resp.train[,j], penalty = "enet", alpha = alpha[i])
      df[j] <- eval(parse(text = paste("fit[k]$", df_result, "$df", sep = "")))
      lambda[j] <- eval(parse(text = paste("fit[k]$", df_result, "$tuning", sep = "")))
      est.err[j] <- (1 / nrow(resp.train)) * 
        sum((predict(fit, obs.train)[, criterion] - resp.train[, j])^2)
      pred.err[j] <- (1 / nrow(resp.val)) * 
        sum((predict(fit, obs.val)[, criterion] - resp.val[, j])^2)
      corr[j] <- cor(predict(fit, obs.val)[, criterion], resp.val[, j])
    }
    ave.df[i] <- mean(as.numeric(df))
    ave.lambda[i] <- mean(lambda)
    ave.est.err[i] <- mean(est.err)
    ave.pred.err[i] <- mean(pred.err)
    ave.corr[i] <- mean(corr)
    lambda <- NULL
    df <- NULL
    est.err <- NULL
    pred.err <- NULL
    corr <- NULL
  }
  result <- data.frame(ave.lambda, ave.est.err, ave.pred.err, ave.corr, ave.df)
  result <- round(t(result), 2)
  colnames(result) <- c("Ridge", "alpha.25", "alpha.5", "alpha.75",
                        "Lasso")
  rownames(result) <- c("Lambda", "Estimation Error", "Prediction Error", 
                        "Correlation Score", "Non-zero coefficient")
  return(result)
}

predValidation <- function(fit, lambda){
  # compute correlation score on a validation set not used for the 
  # cross-validation (CV) or the estimation stability CV (ESCV).
  # Other results are returned by this function to simplify the
  # summary for the results (that's why Train.o and Train.r are inputs)
  # Args:
  #    - Train.o : data frame with the observation training set used for the CV
  #    - Train.r : data frame with the response training set used for the CV
  #    - lambda : numeric lambda to use to fit the model
  # 
  # Output:
  #    - vector with the input lambda, the correlation score, the prediction
  #    error on the validation set and the number of non null coefficient with
  #    lambda.
  corr <- cor(predict(fit, obs.val, s = lambda)[, ,1], resp.val)
  corr <- mean(diag(corr))
  pred.err <- mean(apply((predict(fit, obs.val, s = lambda)[, ,1] - 
                            resp.val)^2, 2, sum))
  return(c(Lambda = lambda, Correlation = corr, 
           Prediction = pred.err, DF = fit$df))  
}


# take correlation instead of mse


cvESCV <- function(lambda, V, alpha){
  # return the results table for cross-validation (CV) and estimation stability 
  # cross-validation (ESCV)
  #
  # Args:
  #    - lambda : vector with the lambda we want to test. This vector has a fixed
  #     sequence allowing us to align all our solutions for each fold.
  #     It indexes the solution path
  #    - V : integer for the number of folds in CV
  #    - alpha : vector of the elastic net mixing parameter we want to test
  #
  # Output:
  #    - data frame with the summary of the results from both CV and ESCV.
  #    The data frame contains lambda, correlation score, prediction and
  #    estimation error and the number of non-null coefficients.
  result <- rep(NA, 5)
  for (j in 1:length(alpha)){
    # alpha varies from ridge to lasso regularization
    # alpha = 0 is ridge, alpha = 1 is lasso
    # 0 < alpha < 1 is elastic net regularization
    my.alpha <- alpha[j]
    
    # build the folds
    N <- nrow(obs.train)
    sub <- list()
    # size is the number of observations in each fold
    size <- trunc(N/V)
    R <- 1:N
    for(i in 1:V){
      sub[[i]] <- sample(R, size)
      R <- R[! R%in% sub[[i]]]
    }
    
    # y.hat is the predicted response. It has 4 dimensions.
    # First dimension is the pictures of a fold, second dimension is the voxels,
    # third dimension is the lambda considered, fourth dimension is the fold considered
    # for the cross-validation
    y.hat <- array(dim = c(size, ncol(resp.train), length(lambda), V))
    
    # mse.hat is mean square error averaged on the voxels with the fold of the cross
    # validation for the rows and the lambda for the columns
    cor.hat <- array(dim = c(V, length(lambda)))
    
    for(k in 1:V){
      # cross-validation on the subset k
      # k th subset removed, fitting on the remaining training set, 
      # prediction on the testing set (ie k th subset)
      Train.o <- obs.train[-sub[[k]], ]
      Test.o  <- obs.train[ sub[[k]], ]
      Train.r <- resp.train[-sub[[k]], ]
      Test.r <- resp.train[ sub[[k]], ]
      Test.r <- array(data = Test.r, dim = c(nrow(Test.r), ncol(Test.r), 
                                             length(lambda)))
      fit <- glmnet(Train.o, Train.r, family = "mgaussian", alpha = my.alpha,
                    lambda = lambda)
      # y.hat has 4 dimensions
      y.hat[, , ,k] <- predict(fit, newx = Test.o)
      # mse is average on the voxels
      #mse.hat[k, ] <- (1/ncol(resp.train)) * colSums(
      #  apply((predict(fit, newx = Test.o) - Test.r)^2, c(2,3), mean))
      for (l in 1:length(lambda)){
        cor.hat[k, l] <- mean(diag(cor(predict(fit, newx = Test.o)[, , l],
                                       Test.r[, , l])))
      }
    }
    # lambda for CV is the lambda for the minimum mse
    cor.hat.ave <- apply(cor.hat, 2, mean)
    lambda.cv <- lambda[which.max(cor.hat.ave)]
    fit.cv <- glmnet(obs.train, resp.train, lambda = lambda.cv, 
                     family = "mgaussian")
    # result cv takes the lambda minimum for cv and use the validation set
    # not used before to compute the correlation
    result.cv <- c(predValidation(fit.cv, lambda.cv), Estimation = max(cor.hat.ave))
    # lambda for ESCV is the lambda for the minimum var(Yhat)/mean(Yhat)
    # average Yhat for lambda fixed
    y.hat.ave.lambda <- apply(y.hat, c(1, 2, 3), mean)
    # norm 2 of the average of Yhat
    y.hat.ave.norm2 <- apply(apply(y.hat.ave.lambda^2, c(2, 3), sum), 2, mean)
    # Yhat average in of dimensions 4
    y.hat.ave.lambda <- array(data = y.hat.ave.lambda, 
                           dim = c(size, ncol(Test.r), length(lambda), V))
    # variance of Yhat as defined in the paper ESCV
    var.y <- apply(apply(apply((y.hat - y.hat.ave.lambda)^2, c(2, 3, 4), sum), 
                         c(1, 2), mean), 2, mean)
    # lambda with the minimum var(Yhat)/mean(Yhat)
    es <- var.y/y.hat.ave.norm2
    lambda.escv <- lambda[which.min(es)]
    # lambda chosen for the ESCV is the minimum of the lambda between the lambda
    # from CV and the lambda for ESCV
    lambda.max <- max(lambda.cv, lambda.escv)
    fit.escv <- glmnet(obs.train, resp.train, lambda = lambda.escv, 
                       family = "mgaussian")
    # result escv takes the lambda minimum for escv and use the validation set
    # not used before to compute the correlation
    result.escv <- c(predValidation(fit.escv, lambda.max), Estimation = min(es))
    # summary of the results
    result <- cbind(result, data.frame(CV = result.cv, ESCV = result.escv))
  }
  return(result[ , 2:(length(alpha) * 2 + 1)])
}



stabModels <- function(fit){
  # compute correlation accross different bootstrap sample to assess he stability
  # of the models
  # 
  # Args:
  #      - fit : output of the function glmnet (object of class glmnet)
  # 
  # Output:
  #      - stability.model : a vector with the correlation for each bootsrap
  #      bootstrap done in parallel
  
  ncores <- 2
  registerDoParallel(ncores)
  # initialize the vector with the final correlation
  corr.boot <- NULL
  stability.model <- foreach(i = 1:100) %dopar% {
    # for each bootstrap sample compute the correlation between the bootstrapped
    # fitted observation and the bootstrapped response (same bootstrap sample)
    boot.index <- sample(1:nrow(resp_dat), nrow(resp_dat), replace = TRUE)
    obs.boot <- obs[boot.index, ]
    resp.boot <- resp[boot.index, 1]
    obs.boot <- cbind(rep(1, nrow(obs)), obs.boot)
    corr.boot[i] <- cor(obs.boot %*% as.numeric(coef(fit)), resp.boot)
    return(corr.boot[i])
  }
  stability.model <- unlist(stability.model)
  return(stability.model)
}


stabFeatures <- function(alpha, voxel){
  # Compute the coefficients accross 1000 bootstrap samples and count the 
  # number of times coefficients were non null
  # 
  # Args:
  #    - alpha : elastic parameter (numeric) for the GLM with regularization
  #    - voxel : integer designing the voxel chosen
  #
  # Output:
  #    - data.frame with the number of times over 1000 bootstraps the coefficients
  #    were non null
  ncores <- 8
  registerDoParallel(ncores)
  coef.boot <- list()
  # bootstrap
  stability.features <- foreach(i = 1:1000) %dopar% {
    compteur <- rep(0, ncol(obs))
    boot.index <- sample(1:nrow(resp_dat), nrow(resp_dat), replace = TRUE)
    obs.boot <- obs[boot.index, ]
    resp.boot <- resp[boot.index, voxel]
    # fit the bootstrap sample
    fit <- cv.glmnet(obs.boot, resp.boot, family = "gaussian", nlambda = 100, 
                     parallel = TRUE, type.measure = "mse", nfolds = 10, 
                     alpha = alpha)
    # find non null coefficient
    coef <- coef(fit)[, 1][2:length(coef(fit)[, 1])]
    index.non.null <- which(coef != 0)
    names(index.non.null) <- NULL
    # create a vector and add 1 each time a feature has a non null coefficient
    compteur[index.non.null] <- 1
    return(compteur)
  }
  return(stability.features) 
}