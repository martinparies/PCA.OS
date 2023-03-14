#' choice.component
#'
#' Helps the user to choose the appropriate number of PCAOS component to optimize
#'
#' @param data a data frame with n rows (individuals) and p columns (numeric, nominal and/or ordinal variables)
#'
#' @param level.scale vector(length p) giving the nature of each variable. Possible values: "nom", "ord", "num"
#'
#' @param nb.comp.to.investigate Number of components to investigate (default=ncol(data))
#'
#' @param supp.var a vector indicating the indexes of the supplementary variables
#'
#' @return Data frame with global Loss value and percentage of the quantified variables inertia explained, for different model dimensionality (H). The last column provides the variation in percentage between to two successive size of dimensionality.
#'
#' @author
#' \itemize{
#'   \item Martin PARIES (Maintainer: \email{martin.paries@oniris-nantes.fr})
#'   \item Evelyne Vigneau
#'   \item Stephanie Bougeard
#' }
#'
#' @examples
#' data("antibiotic")
#'
#' #Construction of the "level.scale" argument for this dataset
#' level.scale <- rep(NA,ncol(antibiotic))
#' level.scale[c(2,3,4)] <- "num"
#' level.scale[c(1,5,6,7,8,9,10,11,12,13,14,15)] <- "nom"
#' level.scale[c(1,15)] <- "ord"
#'
#' res.choice <- choice.component(antibiotic,level.scale)
#' res.choice
#'
#' @export
#'
choice.component <-
  function(data,
           level.scale = rep("num", ncol(data)),
           nb.comp.to.investigate = 5,
           supp.var = NULL) {
    res <- percentage <- NULL
  for (i in 1:nb.comp.to.investigate){
    res.PCAOS <- PCAOS(data,level.scale,nb.comp = i,print = FALSE,init = 'rdm',threshold = 10e-6,supp.var = supp.var)
    res[i] <- res.PCAOS$Algo$loss.tot
    percentage[i] <- res.PCAOS$Dimension.reduction$inertia[i,2]
  }

  pourcentage = rep(NA,nb.comp.to.investigate)
  for (H in 1:nb.comp.to.investigate){
    pourcentage[H] = H * (1 - res[H]) * 100
    delta = rep(NA,nb.comp.to.investigate)
    delta[1] = pourcentage[1]
    for (H in 2:nb.comp.to.investigate){
      delta[H] = pourcentage[H] - pourcentage[H-1]
    }
  }
  data <- data.frame(Component = 1:nb.comp.to.investigate,Loss = round(res,4),Percentage = round(pourcentage,4) ,Improvement = round(delta,4))
  colnames(data) = c("nb.comp","Loss","Percentage (%)","Improvement  (%)")
  choice.inertia <- ggplot2::ggplot(data=data, ggplot2::aes(x= data[,1], y=data[,4])) +
    ggplot2::geom_bar(stat="identity", width=0.5) + ggplot2::theme_classic(base_size = 10) +
    ggplot2::xlab("Number of components in the model") +
    ggplot2::ylab("Improvement of inertia") +
    ggplot2::geom_text(ggplot2::aes(label=paste(round(delta,2),"%")), vjust= -0.1, hjust = 0.3, color="black", size=3.5)
  print(choice.inertia)
  return(data)
}

