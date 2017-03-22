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