#' print.PCAOS
#'
#' Print results of PCAOS method
#'
#' @param x an object of class PCAOS
#'
#' @param file A connection, or a character string naming the file to print to. If NULL (the default), the results are not printed in a file
#'
#' @param sep character string to insert between the objects to print
#' @param ...	further arguments passed to or from other methods
#'
#' @examples
#' data (antibiotic)
#' # Level of scaling of each variable
#' # Manually
#' level.scale <- rep(NA,ncol(antibiotic)) #Setting level.scale argument
#' level.scale[c(3,4)] <- "num"
#' level.scale[c(6:14)] <- "nom"
#' level.scale[c(1,15)] <- "ord"

#' # Or using nature.variables()
#' level.scale <- rep(NA,ncol(antibiotic))
#' res.nature <- nature.variables(antibiotic)
#' level.scale [res.nature$p.numeric] <- "num"
#' level.scale [res.nature$p.quali] <- "nom"
#' #Warning; the ordinal nature of variables can not be detected automaticaly.
#' level.scale[c(1,15)] <- "ord"
#'
#' # PCAOS
#' res.PCAOS <- PCAOS(
#' data = antibiotic,
#' level.scale = level.scale,
#' nb.comp = 2)
#'
#' print(res.PCAOS)
#'
#' @export print.PCAOS
#' @export

print.PCAOS <- function (x, file = NULL, sep = ";", ...){
  res.PCAOS <- x
  if (!inherits(res.PCAOS, "PCAOS")) stop("non convenient data")
  cat('##############################################################')
  cat("\n")
  cat('** The folowing results are available from PCAOS analysis **')

  res <- matrix(NA,14,1)
  res[1, ] <- '$Dimension.reduction$weights: weights of variables'
  res[2, ] <- '$Dimension.reduction$components: components of the PCAOS analysis'
  res[3, ] <- '$Dimension.reduction$inertia : percentage of inertia explained'
  res[4, ] <- '$Quantifications$quantified.data: quantified data trough Optimal Scaling'
  res[5, ] <- '$Quantifications$quantified.categories.nom (resp.ord): quantified categories trough Optimal Scaling'
  res[6, ] <- '$Quantifications$level.scale: levels of scaling or variables'
  res[7, ] <- '$Quantifications$data: original dataset'
  res[8, ] <- '$Algorithm$summary: summary of level of scaling'
  res[9, ] <- '$Algorithm$loss.tot: global loss'
  res[10, ] <- '$Algorithm$stockiter: evolution of the criterion'
  res[11, ] <- '$Supp.var$var.supp: original supplementary variables'
  res[12, ] <- '$Supp.var$level.scale.supp: level of scaling of supplementary variables'
  res[13, ] <- '$Supp.var$coord.supp.num: cordinates of supplementary numeric variables'
  res[14, ] <- '$Supp.var$coord.supp.quali: coordinates of qualitatve variables '
  rownames(res) <- rep('*',14)
  colnames(res) <- ''
  print(res)
  cat('##############################################################')

}
