#' choice.component
#'
#' Helps user choose appropriate number of component to optimise for PCAOS method.
#'
#' @param data a data frame with n rows (individuals) and p columns (numeric, nominal and/or ordinal variables)
#'
#' @param nature vector(length p) giving the nature of each variable. Possible values: "nom", "ord", "num"
#'
#' @param nb.comp.to.investigate Number of components to investigate (default=ncol(data))
#'
#' @return Data frame with global Loss value and percentage of the quantified variables inertia explained, for different model dimensionality (H). The last column provides the variation in percentage between to two successive size of dimensionality.
#'
#' @export
#'
choice.component <- function(data,nature = rep("num",ncol(data)),nb.comp.to.investigate = ncol(data)){
  res <- sapply(1:nb.comp.to.investigate,function(x){PCAOS(data,nature,nb.comp = x)$loss.tot})
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
    ggplot2::geom_bar(stat="identity", width=0.5) + ggplot2::theme_classic(base_size = 15) +
    ggplot2::xlab("Number of components in the model") +
    ggplot2::ylab("Improvement of inertia") +
    ggplot2::geom_text(ggplot2::aes(label=paste(round(delta,2),"%")), vjust= -0.1, color="black", size=5)
  print(choice.inertia)
  return(data)
}

