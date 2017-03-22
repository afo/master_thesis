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