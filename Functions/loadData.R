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
