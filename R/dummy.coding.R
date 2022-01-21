# dummy.coding
#
# dummy.coding
#
# @param data raw variables (vector or matrix)
#
# @return dummy coded variables
#

dummy.coding <- function(data){
  resDIS <- NULL
  modal <- NULL
  for (i in 1:ncol(data)){
    col <- data[,i]
    var <- colnames(data)
    nvar <- var[i]
    nindiv <- length(col)
    col <- as.factor(col)
    resDISvar <- matrix(0, nindiv, length(levels(col)))
    resDISvar[(1:nindiv) + nindiv * (unclass(col) - 1)] <- 1
    rownames(resDISvar) <- rownames(data)
    colnames(resDISvar) <- paste(nvar,levels(col),sep ="/")
    resDIS <- cbind(resDIS,resDISvar)
    modal <- c(modal,levels(col))
  }
  res <- list(data = resDIS,modal = modal)
  return(res)
}
