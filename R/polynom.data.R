# polynom.data
#
# polynom.data
#
# @param var variable to pretreat
# @param D Degree of the relation between quantified variables and components
#
# @return dummy coded variables
#
polynom.data <- function(var,D){
  var[which(is.na(var)),] <- 0
  nbINDIV <- nrow(var)
  nbCOL <- ncol(var)
  nbCOL <- nbCOL * (D+1)
  res <- matrix(0,nbINDIV,nbCOL)

  if (D==1){
    sequencePUISSANCE <- rep(c(0,1),times = nbCOL)
    sequenceRAWdata <- rep(c(1:nbCOL),each=D+1)
    for (i in 1:length(res[1,])){
      res[,i] <- var[,sequenceRAWdata[i]] ^ sequencePUISSANCE[i]
    }
  }

  if (D==2){
    sequencePUISSANCE <- rep(c(0,1,2),times = nbCOL)
    sequenceRAWdata <- rep(c(1:nbCOL),each=D+1)
    for (i in 1:length(res[1,])){
      res[,i] <- var[,sequenceRAWdata[i]] ^ sequencePUISSANCE[i]
    }
  }

  return(res)
}
