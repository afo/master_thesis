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