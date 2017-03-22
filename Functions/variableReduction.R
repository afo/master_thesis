###### VARIABLE REDUCTION & FINAL ALGORITHM ######

#------ Main Functions ------#

getGoodnessOfFitMatrix <- function(Y, traits, vec, n, a){
  # Computes goodness of fit variates (RMSE, AIC and BIC) for k=1:24 where k
  # is the number of regression variables.
  # The values of RMSE, AIC and BIC that are returned is the mean of these 
  # variates over each of the seven linear models given by each of the latent 
  # traits regressed over the item response matrix
  #
  # Args:
  #   Y:        item response matrix
  #   traits:   first and second order trait estimates (from MCMC)
  #   vec:      feature ranking vector (boruta or GBM)
  #   n:        number of cross-validations
  #   a:        percentage of data that is training data in the 
  #             cross-validations
  #
  # Returns:
  #   results:  matrix with the mean RMSE, AIC and BIC-values for each
  #             of k=1:24 regression variables
  #
  
  results <- matrix(nrow=24, ncol=3)
  colnames(results) <- c("RMSE", "AIC", "BIC")
  row.names(results) <- 1:24
  
  for(i in 1:24){
    Y.red <- as.matrix(Y[,vec[1:i]])
    lm.mods <- lapply(1:7, function(j)lmCrossValidationRun(Y.red, traits[,j], n, a))
    mean.RMSE <- mean(sapply(1:7, function(j)lm.mods[[j]]$RMSE))
    mean.AIC <- mean(sapply(1:7, function(j)AIC(lm.mods[[j]]$mod)))
    mean.BIC <- mean(sapply(1:7, function(j)BIC(lm.mods[[j]]$mod)))
    results[i,] <- c(mean.RMSE, mean.AIC, mean.BIC)
  }
  return(results)
}


lmCrossValidationRun <- function(X, Y, n, a){
  # Performs a linear regression with cross validation.
  # 
  # Args:
  #   X:    explanatory variables
  #   Y:    dependent variable
  #   n:    number of cross validations
  #   a:    data split, percentage of data used for training
  # 
  # Returns:
  #   mod:  model with lowest RMSE obtained in CV
  #   RMSE: RMSE-value of best model
  
  mods <- (lapply(1:n, function(i)lmRun(X,Y,a)))
  best.id <- which.min(sapply(1:n, function(i)mods[[i]]$rmse))
  best.mod <- mods[[best.id]]$mod
  RMSE <- mods[[best.id]]$rmse
  
  return(list(mod = best.mod, RMSE = RMSE))
}

getYhat <- function(c, x){
  # Computes estimates of dependent variable from the regression parameters
  # obtained in a linear regression
  # 
  # Args:
  #   c:    vector of regression coefficients
  #   x:    data matrix or row vector (explanatory variables), [n x k] 
  #         where n=samples, k=regression parameters
  #
  # Returns:
  #   yhat: estimates of the dependent variable
  
  n <- dim(as.matrix(x))[1]
  x.full <- cbind(rep(1,n), x)
  yhat <- as.matrix(x.full)%*%as.matrix(c)
  return(yhat)
}

getRMSE <- function(yhat, y){
  # Computes the Root Mean Square Error (RMSE) of a parameter estimate
  #
  # Args:
  #   yhat: parameter estimate
  #   y:    true parameter value
  #
  # Returns:
  #   RMSE: Root Mean Square Error of parameter estimate 'yhat'
  
  y <- as.matrix(y)
  yhat <- as.matrix(yhat)
  RMSE <- sqrt(mean((yhat-y)^2))
  return(RMSE)
}

oneToTenScale <- function(data, lo, hi) {
  # Casts the index over the interval I[1,10]
  #
  # Args:
  #   data:         vector containg data to be cast on new interval
  #   lo:           former lowest possible value, recast as lo -> 1
  #   hi:           former highest possible value, recast as hi -> 10
  #
  # Returns:
  #   recast.data:  data recast on interval I[1,10]
  
  old.range = (hi-lo)  
  new.range = (10-1) 
  
  recast.data <- sapply(1:length(data), function(i)(((data[i]-lo)*new.range)/old.range)+1)
  recast.data[which(recast.data>10)] <- 10
  
  return(recast.data)
}

#------ Auxiliary Functions ------#

lmRun <- function(X, Y, a){
  # Auxiliary function
  #
  # Performs one linear regression with the data set
  # randomly split in training set and test set.
  #
  # Args:
  #   X:    explanatory variables
  #   Y:    dependent variable
  #   a:    data split, percentage of data used for training
  # 
  # Returns:
  #   mod:  linear regression model
  #   rmse: RMSE-value for mod
  
  sets <- getSets(X,Y,a)
  X.train <- sets$x.train
  X.test <- sets$x.test
  Y.train <- sets$y.train
  Y.test <- sets$y.test
  
  lm.mod <- lm(Y.train ~ X.train)
  lm.coef <- lm.mod$coefficients
  yhat <- getYhat(lm.coef, X.test)
  lm.rmse <- getRMSE(yhat, Y.test)
  return(list(rmse = lm.rmse, mod = lm.mod))
}

getSets <- function(X, Y, a){
  # Auxiliary function
  #
  # Split data set in training and test set.
  #
  # Args:
  #   X:    explanatory variables
  #   Y:    dependent variable
  #   a:    data split, percentage of data used for training
  # 
  # Returns:
  #   x.train:  explanatory variables, training set
  #   y.train:  dependent variable, training set
  #   x.test:   explanatory variables, test set
  #   y.test:   dependent variable, test set
  
  n <- dim(X)[1]
  train.rows <- sample(1:n, a*n)
  x.train <- X[train.rows, ]
  x.test <- X[-train.rows, ]
  
  y.train <- Y[train.rows]
  y.test <- Y[-train.rows]
  return(list(x.train = x.train, x.test = x.test, 
              y.train = y.train, y.test = y.test))
}