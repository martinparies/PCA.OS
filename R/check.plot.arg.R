#' check.plot.ar
#'
#' check if data and arguments are good to be plot
#'
#' @param choice the graph to plot possible values are "screeplot","quantif","indiv","cor","modalities","mixed","squared loadings". See Details.
#' @param nature vector(length p) giving the nature of each variable. Possible values: "nom", "ord", "num"
#' @param nature.supp nature of supplementary variables
#' @param supp.var if the user asked to plot the supplementary variable or not
#' @param comp axes to use in the plots
#' @param nb.comp number of components in the PCAOS model
#'
#' @return stop function if argument are wrong
#'
check.plot.arg <- function(choice,nature,nature.supp,supp.var, comp, nb.comp){
  if (!(choice %in% c("screeplot","numeric","qualitative","quantif","ind","quali.supp","mixed"))){
    stop("Values of choice should be one of : screeplot,numeric,qualitative,quantif,ind,quali.supp,mixed")
  }
  if (choice == "qualitative"){
    if (!any(nature == "nom")){
      stop("No variable defined as qualitative (i.e nominal or ordinal")
    }
    if(supp.var == TRUE & !any(nature.supp == "nom" | nature.supp == "ord")){
      stop("Supplementary variable is not qualitative")
    }
  }

  if (choice == "numeric"){
    if (!any(nature == "nom")){
      stop("No variable defined as numeric")
    }
    if(supp.var == TRUE & !any(nature.supp == "num")){
      stop("Supplementary variable is not numeric")
    }
  }

  if (choice == "mixed"){
    if ((!any(nature == "nom") | !(any(nature == "ord"))) & !any(nature == "num") ){
      stop("Variables are not mixed")
    }
  }

  if (any(comp > nb.comp)){
    stop("Wrong axes to plot")
  }

}
