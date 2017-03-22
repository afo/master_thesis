###### ALL FUNCTIONS ######


###### LOAD & PREPARE REAL DATA ######

#------ Main Function ------#

readData = function(reverse = TRUE, NAThree = TRUE, clean = TRUE, onlyQT = TRUE) {
  # Loads and returns item response data set
  #
  # Args: 
  #   reverse: all items have best possible response k=5
  #   NAThree: NA-values replaced with k=3
  #   clean:   data set cleaned from bad samples
  #   onlyQT:  only return QT-items, i.e. not CZ, CF or SDR
  #
  # Returns:
  #   data:    data set as specified by input
  
  if (clean) {data.full <- returnCleanData()}
  else {data.full = read.csv("dat.csv")}
  
  
  if (NAThree) {
    data.full[is.na(data.full)] <- 3 
  }
  
  if (reverse) {
    vec = c("QT2","QT4","QP1","QP3")
    
    for (i in 1:length(vec)) {
      tmp <- data.full[,vec[i]]
      tmp = -(tmp - 3)+3
      data.full[,vec[i]] = tmp
    }
  }
  
  if (onlyQT){
    data.full = data.full[,8:31]
  }
  else {
    data.full = cbind(data.full[,2:4], data.full[,8:31])
  }
  
  return(data.full)
}

#------ Auxiliary Functions ------#

returnCleanData <- function(){
  # Auxiliary function
  #
  # Load data from file "dat.csv" and clean it from bad samples.
  # 
  # Returns:
  #   data: reduced data set with only item responses, [878x24]
  
  dfull = read.csv("dat.csv")
  
  data<-cleanData(dfull)
  
  start<-grep("QT1",names(data))
  end <- grep("QC4",names(data))
  answers <- data[,start:end]
  
  for (i in 1:5) {
    repetitions<-rowMatch(answers,rep(i,24))
    if (!is.null(repetitions)) {
      data <- data[-repetitions,]
    }
  }
  
  return(data)
}

rowMatch <- function(df, vec){
  # Auxiliary function
  #
  # Identify item response vectors identical to 'vec'
  #
  # Args:
  #   Y: item response matrix
  #
  # Returns:
  #   ind: indices for samples where repsonses to all items are
  #        equal
  
  vec = as.integer(vec)
  rows <- dim(df)[1]
  ind <- NULL
  for (i in 1:rows) {
    if (isTRUE(all.equal(as.integer(df[i,]),vec))) {
      ind = c(ind,i)
    }
  }
  return(ind)
}

cleanData <- function(df) {
  # Auxiliary function
  #
  # Identify and remove all but first sample submitted by the same individual
  #
  # Args:
  #   df: full data set
  #
  # Returns:
  #   df: full data set cleaned from duplicate samples
  
  # Identify duplicates
  mails <- df[,7]
  dups <- mails[duplicated(mails)]
  dups <- dups[-which(dups=="")]
  dups <- unique(dups)
  
  # Remove all but first occurance from data set
  n = length(dups)
  for (i in 1:n) {
    ind<-grep(dups[i],mails)
    mails <- mails[-ind[-1]]
    df <- df[-ind[-1],]
  }
  return(df)
}


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


###### POSTERIOR PREDICTIVE CHECK ######

#------ Main Functions ------#

meanPPCM <- function(jagsUI.model, Y){
  # Runs a Posterior Predictive Mean Check
  #
  # Args:
  #   jagsUI.model: Model created from JagsUI-package with replicated data
  #   Y: Item response matrix
  #
  # Returns:
  #   p/pLo/pHi: mean p-values for the PPC
  
  dims <- dim(jagsUI.model$sims.list$theta)
  
  data.means <- matrix(rep(rowMeans(Y),dims[1]), nrow = dims[1], ncol = dims[2], byrow=TRUE)
  sim.means <- t(apply(jagsUI.model$sims.list$Ynew,1,rowMeans))
  
  
  hist(sim.means,main="",xlab="")
  abline(v=mean(rowMeans(Y)),lwd=4)
  title(xlab=expression(paste("Row Means, ", 1/n, sum(x[i]^new, i==1, n))), line=3.5)
  
  p <- mean(sim.means==data.means)
  pLo <- mean(sim.means < data.means)
  pHi <- mean(sim.means > data.means)
  
  return(c(p,pLo,pHi))
}

stdPPCM <- function(jagsUI.model, Y){
  # Runs a Posterior Predictive Std Check
  #
  # Args:
  #   jagsUI.model: Model created from JagsUI-package with replicated data
  #   Y: Item response matrix
  #
  # Returns:
  #   p/pLo/pHi: std p-values for the PPC
  
  dims <- dim(jagsUI.model$sims.list$theta)
  
  data.sds <- matrix(rep(rowSds(as.matrix(Y)),dims[1]), nrow = dims[1], ncol = dims[2], byrow=TRUE)
  sim.sds <- t(apply(jagsUI.model$sims.list$Ynew,1,rowSds))
  
  
  hist(sim.sds,main="",xlab="", xlim=c(0,2))
  abline(v=mean(rowSds(as.matrix(Y))),lwd=4)
  title(xlab=expression(paste("Row Standard Deviations, ", sigma[s](x[si]^rep))), line=3.5)
  
  p <- mean(sim.sds == data.sds)
  pLo <- mean(sim.sds < data.sds)
  pHi <- mean(sim.sds > data.sds)
  
  return(c(p,pLo,pHi))
}


###### PLOT ITEM CHARACTERISTIC CURVES ######

#------ Main Function ------#

plotICC <- function(mod){
  # Plots the HO-GRM Item Characteritic Curves for each item.
  # Only works for data structures equal to the real data.
  # 
  # Args:
  #   mod: jagsUI-model of real data set
  
  a <- mod$mean$a
  b <- mod$mean$b
  theta <- seq(-4,4,0.05)
  
  par(mfrow=c(1,1))
  titles <- c("QT1","QT2", "QT3", "QT4", "QF1", "QF2", "QF3", "QF4", 
              "QD1", "QD2", "QD3", "QD4", "QB1", "QB2", "QB3", "QB4", 
              "QP1", "QP2", "QP3", "QP4", "QC1", "QC2", "QC3", "QC4")
  subtitle <- c("Trust", "Resilience", "Diversity", "Belief", "Perfection", "Collaboration")
  
  count = 0
  
  for(i in 1:length(a)){
    if((i-1) %% 4 == 0){
      count = count + 1
    }
    a_temp <- a[i]
    b_temp <- b[i,]
    probs <- getItemICC(a_temp,b_temp,theta)
    plotProbs(probs, theta, titles[i], subtitle[count], count)
  }
}


#------ Auxiliary Functions ------#

getItemICC <- function(a, b, x){
  # Auxiliary function
  #
  # Get item characteristic curve for each possible item response for one item
  #
  # Args:
  #   a: item discrimination parameter    
  #   b: item difficulty parameter
  #   x: vector with values over the latent trait continuum
  #
  # Returns:
  #   p: item characteristic cruves for all possible item responses 
  
  logit <- function(y){exp(y)/(1+exp(y))}
  
  p1 <- 1-logit(a*(x-b[1]))
  p2 <- logit(a*(x-b[1])) - logit(a*(x-b[2]))
  p3 <- logit(a*(x-b[2])) - logit(a*(x-b[3]))
  p4 <- logit(a*(x-b[3])) - logit(a*(x-b[4]))
  p5 <- logit(a*(x-b[4]))
  
  p <- cbind(p1,p2,p3,p4,p5)
  return(p)
} 

plotProbs <- function(p, x, title, subtitle, count){
  # Auxiliary function
  #
  # Plots the Item Characteritic Curves for one item
  #
  # Args:
  #   p:        item characteristic curves for item
  #   x:        vector with values over the latent trait continuum
  #   title:    item name
  #   subtitle: domain name
  #   count:    domain number
  
  colors <- c("blue", "red", "green", "orange", "black")
  
  for(i in 1:dim(p)[2]){
    if(i==1){
      plot(x, p[,i], type="l", col=colors[i], xlab="" , ylab="",ylim=c(0,1) )
      title(title, line=1.6)
      title(ylab=substitute(paste("P(",x[ik],"|",theta[count]^(1),")", sep=""),
                            list(count=count)), mgp=c(2.5,2,0))
      title(xlab=substitute(paste(subtitle, ", ", theta[count]^(1)),
                            list(subtitle=subtitle, count=count)), line=2.5)
    } else{
      lines(x, p[,i], col=colors[i])
    }
  }
  legend("bottomleft", c("k=1","k=2","k=3","k=4","k=5"), text.width = 0.5, 
         cex=1, col=colors, pch=NA, lwd=1)
}


###### PLOT ITEM INFORMATION CURVES ######

#------ Main Function ------#

plotIIC <- function(mod){
  # Plots the HO-GRM Item Information Curves for each item 
  # as well as Test Information Curves for each domain.
  # Only works for data structures equal to the real data.
  # 
  # Args:
  #   mod: jagsUI-model of real data set 
  
  
  # Parameters
  a <- mod$mean$a
  b <- mod$mean$b
  theta <- seq(-4,4,0.05)
  
  # Plot specifications
  titles = c("Trust", "Resilience", "Diversity", "Belief", "Perfection", "Collaboration")
  subtitles <- c("QT1","QT2", "QT3", "QT4", "QF1", "QF2", "QF3", "QF4", "QD1", "QD2", "QD3", "QD4",
                 "QB1", "QB2", "QB3", "QB4", "QP1", "QP2", "QP3", "QP4", "QC1", "QC2", "QC3", "QC4")
  
  # Get information values for all items over latent trait continuum
  A <- matrix(nrow = length(theta), ncol = 24)
  for(i in 1:length(a)){
    A[,i] <- getItemInfo(a[i], b[i,], theta)
  }
  
  # Plot individual item information
  par(mfrow=c(1,1))
  
  for(q in 1:6){
    j <- (q-1)*4+1
    A.temp <- A[,(j:(j+3))]
    title <- titles[q]
    subtitle <- subtitles[j:(j+3)]
    plotInfo(A.temp,theta,title,subtitle,q)
  }
  
  # Plot test information
  for(q in 1:6){
    j <- (q-1)*4+1
    A.temp <- A[,(j:(j+3))]
    title <- titles[q]
    A.test <- rowSums(A.temp)
    plotTestInfo(A.test, theta, q, titles)
  }
}

#------ Auxiliary Functions ------#

getItemInfo <- function(a, b, x){
  # Auxiliary function
  #
  # Get information values for one item
  #
  # Args:
  #   a: item discrimination parameter    
  #   b: item difficulty parameter
  #   x: vector with values over the latent trait continuum
  #
  # Returns:
  #   A: item information values for each value of vector x 
  
  logit <- function(y){exp(y)/(1+exp(y))}
  
  p1 <- 1
  p2 <- logit(a*(x - b[1]))
  p3 <- logit(a*(x - b[2]))
  p4 <- logit(a*(x - b[3]))
  p5 <- logit(a*(x - b[4]))
  
  A1 <- a^2*(p1*(1-p1)-p2*(1-p2))^2/p1
  A2 <- a^2*(p2*(1-p2)-p3*(1-p3))^2/p2
  A3 <- a^2*(p3*(1-p3)-p4*(1-p4))^2/p3
  A4 <- a^2*(p4*(1-p4)-p5*(1-p5))^2/p4
  A5 <- a^2*(p5*(1-p5))^2/p5
  
  A <- rowSums(cbind(A1,A2,A3,A4,A5))
  return(A)
}

plotInfo <- function(A, x, title, subtitle, q){
  # Auxiliary function
  #
  # Plots the item information for each item in domain q
  #
  # Args:
  #   A:        item information matrix for items in domain q
  #   x:        vector with plotting values (plotting interval)
  #   title:    name of domain q
  #   subtitle: vector with item names for items in domain q
  #   q:        domain
  
  colors <- c("black", "red", "green", "blue")
  
  for(j in 1:4){
    if(j==1){
      title=title
      plot(x,A[,j],type="l", lwd=2, col=colors[j], ylim=c(0,(max(A)+0.2)), main=substitute(bold(paste("Item Information, ", I[i](theta),": ", title))), xlab="", ylab="")
      title(ylab = substitute(paste("Item Information, ", I[i](theta))), mgp=c(2.5,2,0))
      title(xlab = substitute(paste(sub, ", ", theta[count]^(1)),list(sub=title, count=q)), line=2.5)
    } else{
      lines(x,A[,j], type="l", lwd=2, col=colors[j])
    }
  }
  legend("topright", subtitle, cex=1, col=colors, lty=1, text.width = 0.6, lwd=2)
}

plotTestInfo <- function(A, x, q, titles){
  # Auxiliary function
  #
  # Plots the test information for domain q
  #
  # Args:
  #   A:      item information matrix for items in domain q
  #   x:      vector with plotting values (plotting interval)
  #   q:      domain
  #   titles: list of domain names
  
  colors = c("blue", "red", "black")
  
  if(q == 1){
    plot(x,A,type="l", lwd=2, col=colors[1] , ylim=c(0,5), main=expression(bold(paste("Test Information, ", I[j](theta), ": Trust, Resilience, Diversity"))), xlab="", ylab="")
    title(xlab=expression(paste("First order trait, ", theta[j]^(1)), line=0))
    title(ylab=expression(paste("Test Information, ",I[j](theta))), mgp=c(2.5,2,0))
    legend("topright", titles[1:3], cex=1, text.width = 1.4, col=colors, lty=1)
  }
  else if(q==4){
    plot(x,A,type="l", lwd=2, col=colors[1] , ylim=c(0,5), main=expression(bold(paste("Test Information, ", I[j](theta), ": Belief, Perfection, Collaboration"))), xlab="", ylab="")
    title(xlab=expression(paste("First order trait, ", theta[j]^(1)), line=0))
    title(ylab=expression(paste("Test Information, ",I[j](theta))), mgp=c(2.5,2,0))
    legend("topright", titles[4:6], cex=1, text.width = 1.8, col=colors, lty=1)
  }
  else {
    if(q==5 | q==6){
      q=q-3
    }
    lines(x,A,type="l", lwd=2, col=colors[q])
  }
}


###### PLOT CORRELATIONS ######

#------ Main Functions ------#

plotCorVar <- function(Y, theta){
  # Plots the correlation between item response vectors and the
  # second order latent trait together with the variance of each item 
  # response vectors.
  #
  # Args:
  #   Y:      item response vector
  #   theta:  second order latent trait estimates
  
  cor.vec <- cor(Y,theta)
  var.vec <- diag(var(Y))
  
  # Plot
  d = data.frame(x =seq(1,24), n = var.vec, cor.vec = cor.vec)
  par(mar = c(5,5,2,5))
  with(d, plot(x, cor.vec, type="l", lty=1, col="black", 
               ylab=expression(paste(rho(X[i],theta^(2)))),
               ylim=c(0,0.6), xlab=""))
  title(xlab="Item, i", line=2)
  par(new = T)
  with(d, plot(x, n, type = "l", lty=4, axes=F, xlab=NA, ylab=NA, cex=1.2))
  axis(side = 4)
  mtext(side = 4, line = 3, expression(paste(Var(X[i]))))
  legend("bottomleft",
         legend=c("Correlation", "Item variance"),
         lty=c(1,4), col=c("black", "black"), cex=0.9, text.width = 5)
}

plotCorMean <- function(Y, theta){
  # Plots the correlation between item response vectors and the
  # second order latent trait together with the mean item 
  # response for each item.
  #
  # Args:
  #   Y:      item response vector
  #   theta:  second order latent trait estimates
  
  cor.vec <- cor(Y,theta)
  
  # Plot
  d = data.frame(x =seq(1,24), n = colSums(Y)/dim(Y)[1], cor.vec = cor.vec)
  par(mar = c(5,5,2,5))
  with(d, plot(x, cor.vec, type="l", lty=1, col="black", 
               ylab=expression(paste(rho(X[i],theta^(2)))),
               ylim=c(0,0.6), xlab=""))
  title(xlab="Item, i", line=2)
  par(new = T)
  with(d, plot(x, n, type="l", lty=4, pch=1, axes=F, xlab=NA, ylab=NA, cex=1.2))
  axis(side = 4)
  mtext(side = 4, line = 3, expression(paste("Mean Item Response, ", mu[X[i]] )))
  legend("bottomleft",
         legend=c("Correlation", "Mean item resp"),
         lty=c(1,4), col=c("black", "black"), cex=0.9, text.width = 6)
}


###### FEATURE RANKING ######

#------ Main Function ------#

runBoruta <- function(Y, trait) {
  # Runs and plots a boruta run.
  # Replaces all -Inf with 0 in the resulting boruta.
  #
  # Args:
  #   Y:      dependent variable
  #   trait:  explanatory variables
  #
  # Returns:
  #   boruta.train: result from boruta method (including information matrix)
  
  boruta.train <- Boruta(Y, trait, doTrace = 2,maxRuns = 200)
  boruta.train$ImpHistory[which(boruta.train$ImpHistory == -Inf)] <- 0
  
  # Plot
  plot(boruta.train, xlab = "", xaxt = "n")
  lz<-lapply(1:ncol(boruta.train$ImpHistory),function(i)boruta.train$ImpHistory[is.finite(boruta.train$ImpHistory[,i]),i])
  names(lz) <- colnames(boruta.train$ImpHistory)
  Labels <- sort(sapply(lz,median))
  axis(side = 1,las=2,labels = names(Labels), at = 1:ncol(boruta.train$ImpHistory), cex.axis = 0.7)
  
  return(boruta.train)
}


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
