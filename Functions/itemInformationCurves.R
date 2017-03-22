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