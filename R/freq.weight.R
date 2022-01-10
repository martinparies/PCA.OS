#' freq.weight
#'
#' freq.weight
#'
#' @param dummy dummy coded variables
#'
#' @return pretreated dummy coded variables by khi-deux metrics
#'
freq.weight <- function(dummy){
  #standardisation par la proportion
  N <- rep (length(dummy[,1]),length(dummy[1,]))
  effQUALI <- colSums(dummy)
  prop <- effQUALI / N
  dummy <- sweep(dummy,2,prop,"-")
  #Reduction des variables qualitatives par la proportion
  dummy.weighted <- sweep(dummy,2,sqrt(prop),"/")
  return(dummy.weighted)
}
