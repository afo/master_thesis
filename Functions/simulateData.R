###### SIMULATE DATA ######

#------ Main Function ------#

simulateData <- function(n=1000, m=24, d=6, k=5, test=NA, grm=TRUE, seed=71){
  # Simulate a data set according to the HO-GRM or HO-GPCM model structure.
  # It is only possible to simulate data sets where each item has the same
  # number of possible item responses.
  # Returns a data set of simulated item responses as well as the model
  # parameters that defines the model structure.
  # 
  # Args:
  #   n:    number of samples in data set
  #   m:    number of items
  #   d:    number of domains/tests/first-order traits
  #   k:    number of possible item responses (for all items)
  #   test: a vector defining which items belongs to which test, if 
  #         NULL it is assumed the number of items per test is equal
  #         for all items
  #   grm:  if TRUE simulate data according to HO-GRM model structure,
  #         if FALSE simulate HO-GPCM data
  #   seed: random seed number
  #
  # Returns:
  #   Y:       item response matrix, [n x m] 
  #   theta:   second order trait (model parameter), [n x 1]
  #   gamtes:  first order traits (model parameter), [n x d]
  #   epsilon: error term (model parameter), [n x d]
  #   alpha:   item discrimination (model parameter), [m x 1]
  #   beta:    item difficulty (model parameter), [m x (k-1)]
  #   lambda:  factor loading (model parameter), [d x 1]
  
  set.seed(seed)
  
  # Error handling
  if (m %% d != 0 & d != 1 & is.na(test)){
    stop('m/d is not an integer and custom test not defined')
  }
  if (d > m) {stop('More domains than items')}
  
  # Simulate model parameters
  
  theta <- rnorm(n,0,1)
  lambda <- seq(0.6,0.9,(0.9-0.6)/(d-1))
  tau <- (1-lambda^2)
  
  epsilon <- sapply(1:d, function(i)rnorm(n,0,tau[i]))
  gamtes <- sapply(1:d, function(i)lambda[i]*theta+epsilon[,i])
  
  alpha <- rnorm(m,1.8,0.25)
  beta0 <- runif(m,-3,0)
  beta <- cbind(beta0, beta0 + 1, beta0 + 2, beta0 + 3)
  
  # Define which domains the items belong to
  
  if (is.na(test)){
    test <- sort(rep(1:d, m/d))
  }
  
  # Get response matrix
  
  if (grm){
    Y <- getResponseMatrixGRM(alpha, beta, gamtes, test)
  }
  else{
    Y <- getResponseMatrixGPCM(alpha, beta, gamtes, test)
  }
  
  return(list(Y = Y, theta = theta, gamtes = gamtes, epsilon = epsilon, 
              alpha = alpha, beta = beta, lambda = lambda))
}

#------ Auxiliary Functions ------#

getResponseMatrixGRM <- function(alpha, beta, gamtes, test){
  # Auxiliary function
  #
  # Simulate HO-GRM response matrix
  #
  # Args:
  #   alpha:   item discrimination (model parameter), [m x 1]
  #   beta:    item difficulty (model parameter), [m x (k-1)]
  #   gamtes:  first order traits (model parameter), [n x d]
  #   test:    a vector defining which items belongs to which test, if 
  #            NULL it is assumed the number of items per test is equal
  #            for all items
  #
  # Returns: 
  #   Y:      item response matrix, [n x m]
  
  # Initialize parameter values
  n <- dim(gamtes)[1]
  d <- dim(gamtes)[2]
  m <- length(alpha)
  k <- dim(beta)[2]+1
  Y <- matrix(nrow=n, ncol=m)
  
  
  # Get probabilities for each item response k, p[[k]], for each
  # subject, n, and item, m
  logit <- function(t){exp(t)/(1+exp(t))}
  
  t_init <- matrix(ncol=m,nrow=n)
  t <- rep(list(t_init),k)
  t[[1]] <- matrix(rep(1,m*n),ncol=m,nrow=n)
  
  for (i in 2:k){
    for (l in 1:m){
      t[[i]][,l] <- logit(alpha[l]*(gamtes[,test[l]]-beta[l,(i-1)]))
    }
  }
  
  p <- t
  p[[5]] <- t[[5]]
  for (i in 1:(k-1)){
    p[[i]] <- t[[i]] - t[[(i+1)]]
  }
  
  
  # Get response matrix, Y, given cumulative probabilities p 
  
  rand.number <- matrix(runif(n*m), nrow=n, ncol=m)
  Y <- matrix(nrow=n, ncol=m)
  
  cum.prob <- p
  
  for (i in 2:k){
    cum.prob[[i]] <- cum.prob[[i-1]] + p[[i]]
  }
  
  for (i in 1:n){
    for (j in 1:m){
      probs <- rep(0,5)
      for (l in 1:k){
        probs[l] <- cum.prob[[l]][i,j]
      }
      Y[i,j] <- min(which(rand.number[i,j] < probs))
    }
  }
  return(Y)
  
}

getResponseMatrixGPCM <- function(alpha, beta, gamtes, test){
  # Auxiliary function
  #
  # Simulates HO-GPCM response matrix
  #
  # Args:
  #   alpha:   item discrimination (model parameter), [m x 1]
  #   beta:    item difficulty (model parameter), [m x (k-1)]
  #   gamtes:  first order traits (model parameter), [n x d]
  #   test:    a vector defining which items belongs to which test, if 
  #            NULL it is assumed the number of items per test is equal
  #            for all items
  #
  # Returns: 
  #   Y:      item response matrix, [n x m]
  
  # Initialize parameter values
  
  n <- dim(gamtes)[1]
  d <- dim(gamtes)[2]
  m <- length(alpha)
  k <- dim(beta)[2]+1
  Y <- matrix(nrow=n, ncol=m)
  
  
  # Get probabilities for each item response k, p[[k]], for each
  # subject, n, and item, m
  
  t_init <- matrix(ncol=m,nrow=n)
  t <- rep(list(t_init),k)
  t[[1]] <- matrix(rep(1,m*n),ncol=m,nrow=n)
  
  for (l in 1:m) {
    t[[2]][,l]<-exp(alpha[l]*(gamtes[,test[l]]-beta[l,1]))
  }
  
  for (i in 3:length(t)){
    for (l in 1:m) {
      t[[i]][,l] <- t[[i-1]][,l]*exp(alpha[l]*(gamtes[,test[l]]-beta[l,(i-1)]))
    }
  }
  
  sum <- matrix(rep(0,m*n),ncol=m,nrow=n)
  
  for (i in 1:length(t)) {
    sum <- sum + t[[i]]
  }
  
  p <- t
  
  for (i in 1:k) {
    p[[i]] <- t[[i]]/sum
  }
  
  # Get response matrix, Y, given cumulative probabilities p 
  
  rand.number <- matrix(runif(n*m), nrow=n, ncol=m)
  Y <- matrix(nrow=n, ncol=m)
  
  cum.prob <- p
  
  for (i in 2:k) {
    cum.prob[[i]] <- cum.prob[[i-1]] + p[[i]]
  }
  
  for (i in 1:n){
    for (j in 1:m){
      probs <- rep(0,5)
      for (l in 1:k){
        probs[l] <- cum.prob[[l]][i,j]
      }
      Y[i,j] <- min(which(rand.number[i,j] < probs))
    }
  }
  return(Y)
  
}
