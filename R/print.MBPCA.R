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
#'data('antibiotic')
#'antb.uses <- antibiotic[,c('Atb.conso','Atb.Sys')]
#'health <- antibiotic[,c('Age','Loss')]
#'vet.practices <- antibiotic[,c(6:15)]
#'antibiotic <- data.frame(antb.uses,health,vet.practices)


#'# Defining blocks
#'blocks.name =  c("antibiotic.uses","Health.of.turkeys","Veterinary.practices")
#'blocks <- c(2,2,10)
#'
#'# Level of scaling
#'level.scale <- rep(NA,ncol(antibiotic))
#'res.nature <- nature.variables(antibiotic)
#'level.scale [res.nature$p.numeric] <- "num"
#'level.scale [res.nature$p.quali] <- "nom"
#'#Warning; the ordinal nature of variables can not be detected automaticaly.
#'level.scale[c(1,14)] <- "ord"
#'
#' # MBPCAOS
#'res.MBPCAOS <- MBPCAOS(data = antibiotic,
#'                      level.scale = level.scale,
#'                       blocks = blocks,
#'                       blocks.name = blocks.name,
#'                       nb.comp = 3)
#'
#' print(res.MBPCAOS)
#'
#' @export print.MBPCAOS
#' @export

print.MBPCAOS <- function (x, file = NULL, sep = ";", ...){
  res.MBPCAOS <- x
    if (!inherits(res.MBPCAOS, "MBPCAOS")) stop("non convenient data")
  cat('##############################################################')
  cat("\n")
  cat('** The folowing results are available from MBPCAOS analysis **')
  res <- matrix(NA,18,1)
  res[1, ] <- '$Dimension.reduction$weights: weights of variables'
  res[2, ] <- '$Dimension.reduction$components: components of the PCAOS analysis'
  res[3, ] <- '$Dimension.reduction$inertia : percentage of inertia explained'
  res[4, ] <- '$Quantifications$quantified.data: quantified data trough Optimal Scaling'
  res[5, ] <- '$Quantifications$quantified.categories.nom (resp.ord): quantified categories trough Optimal Scaling'
  res[6, ] <- '$Quantifications$level.scale: levels of scaling or variables'
  res[7, ] <- '$Quantifications$data: original dataset'
  res[8, ] <- '$Blocks$block.components: components associated with each block'
  res[9, ] <- '$Blocks$block.weight: weights associated with each block'
  res[10, ] <- '$Blocks$blocks: number of variable in each block'
  res[11, ] <- '$Blocks$blocks.name: name of each block'

  res[12, ] <- '$Algorithm$summary: summary of level of scaling'
  res[13, ] <- '$Algorithm$loss.tot: global loss'
  res[14, ] <- '$Algorithm$stockiter: evolution of the criterion'
  res[15, ] <- '$Supp.var$var.supp: original supplementary variables'
  res[16, ] <- '$Supp.var$level.scale.supp: level of scaling of supplementary variables'
  res[17, ] <- '$Supp.var$coord.supp.num: cordinates of supplementary numeric variables'
  res[18, ] <- '$Supp.var$coord.supp.quali: coordinates of qualitatve variables '

  rownames(res) <- rep('*',18)
  colnames(res) <- ''
  print(res)
  cat('##############################################################')

}
