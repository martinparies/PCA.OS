#' mix.data
#'
#' Splits a mixed data matrix according to the nature of variables. Columns of class integer are considered quantitative. Character of factor are consider qualitative.
#'
#' @param data a data frame with n rows (individuals) and p columns (numeric, nominal and/or ordinal variables)
#'
#' @return
#' \itemize{
#'   \item data.numeric : data.frame with numeric variables
#'   \item p.numeric : columns index of numeric variables
#'   \item nb.numeric : number of numeric variables
#'   \item data.quali : data.frame with numeric variables
#'   \item p.quali : columns index of numeric variables
#'   \item nb.quali :  number of numeric variables
#' }
#'
#' @examples
#' data (antibiotic)
#' res.mix <- PCA.OS::mix.data(antibiotic)
#' res.mix$p.numeric
#' head(res.mix$data.numeric)
#' res.mix$p.quali
#' head(res.mix$data.quali)
#' nature <- rep(NA,ncol(antibiotic)) #Setting nature argument for PCAOS function
#' nature[res.mix$p.numeric] <- "num"
#' nature[res.mix$p.quali] <- "nom"
#' nature[c(1,15)] <- "ord"
#' nature
#' @export
#'
mix.data <- function(data){
  if(!is.data.frame(data)){
    stop("Data should be a data.frame")
  }
  class.data <- lapply(data,class)
  class.data <- unlist(class.data)
  quanti <- which(class.data %in% c("numeric", "integer"))
  quali <- which(class.data %in% c("factor", "character"))
  data.num <- data[,quanti,drop = FALSE]
  data.quali <- data[,quali,drop = FALSE]
  nb.numeric <- length(quanti)
  nb.quali <- length(quali)
  res <-
    list(
      data.numeric = data.num,
      p.numeric = quanti,
      nb.numeric = nb.numeric,
      data.quali = data.quali,
      p.quali = quali,
      nb.quali = nb.quali
    )
  return(res)
}
