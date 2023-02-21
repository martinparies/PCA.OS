#' choice.component
#'
#' Helps the user to choose the appropriate number of MBPCAOS component to optimize
#'
#' @param data a data frame with n rows (individuals) and p columns (numeric, nominal and/or ordinal variables)
#'
#' @param level.scale vector(length p) giving the nature of each variable. Possible values: "nom", "ord", "num"
#'
#' @param blocks vector(length k) with number of variables in each bloc
#'
#' @param blocks.name vector(length k) with names of each bloc
#'
#' @param block.scaling scaling applied to each block. Possible value are : \itemize{
#'   \item "inertia"(default): each quantified block is divided by its total inertia (sum of square).
#'   \item "lambda1" : each quantified block is divided by its the first singular value.
#'   \item "null" : no scaling is applied
#' }
#'
#' @param nb.comp.to.investigate Number of components to investigate (default=ncol(data))
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
#' data('antibiotic')
#' antb.uses <- antibiotic[,c('Atb.conso','Atb.Sys')]
#' health <- antibiotic[,c('Age','Loss')]
#' vet.practices <- antibiotic[,c(6:15)]
#' antibiotic.MB <- data.frame(antb.uses,health,vet.practices)
#'
#'# Defining the blocks
#' blocks.name =  c("antibiotic.uses","Health.of.turkeys","Veterinary.practices")
#' blocks <- c(2,2,10)
#'
#'# Level of scaling
#' level.scale.MB <- rep(NA,ncol(antibiotic.MB))
#' res.nature <- nature.variables(antibiotic.MB)
#' level.scale.MB [res.nature$p.numeric] <- "num"
#' level.scale.MB [res.nature$p.quali] <- "nom"
#' #Warning; the ordinal nature of variables can not be detected automaticaly.
#' level.scale.MB[c(1,14)] <- "ord"
#'
#'# Choice of number of components
#' help(choice.component.MB)
#' res.choice.MB <- choice.component.MB(antibiotic.MB,
#'                                      level.scale.MB,
#'                                      blocks,
#'                                      blocks.name,
#'                                      block.scaling = 'inertia')
#'res.choice.MB
#'
#' @export
#'
#'

choice.component.MB <-
  function(data,
           level.scale,
           blocks,
           blocks.name,
           nb.comp.to.investigate = 5,
           block.scaling = 'inertia') {
    percentage <-
      sapply(1:nb.comp.to.investigate, function(x) {
        PCA.OS::MBPCAOS(data = data,
                level.scale = level.scale,
                blocks = blocks,
                blocks.name = blocks.name,
                nb.comp = x,
                print = FALSE,
                maxiter = 50,
                threshold = 10e-5,
                block.scaling = block.scaling)$Dimension.reduction$inertia[x,2]
      })

    pourcentage = rep(NA,nb.comp.to.investigate)

    for (H in 1:nb.comp.to.investigate){
      pourcentage[H] = percentage[H]
      delta = rep(NA,nb.comp.to.investigate)
      delta[1] = pourcentage[1]
      for (H in 2:nb.comp.to.investigate){
        delta[H] = pourcentage[H] - pourcentage[H-1]
      }
    }

    data <- data.frame(Component = 1:nb.comp.to.investigate,Percentage = round(pourcentage,4) ,Improvement = round(delta,4))
    colnames(data) = c("nb.comp","Percentage (%)","Improvement  (%)")
    choice.inertia <- ggplot2::ggplot(data=data, ggplot2::aes(x= data[,1], y=data[,3])) +
      ggplot2::geom_bar(stat="identity", width=0.5) + ggplot2::theme_classic(base_size = 10) +
      ggplot2::xlab("Number of components in the model") +
      ggplot2::ylab("Improvement of inertia") +
      ggplot2::geom_text(ggplot2::aes(label=paste(round(delta,2),"%")), vjust= -0.1, hjust = 0.3, color="black", size=3.5)
    print(choice.inertia)
    return(data)
  }
