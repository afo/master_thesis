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
