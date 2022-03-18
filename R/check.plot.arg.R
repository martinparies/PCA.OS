# check.plot.ar
#
# check if data and arguments are good to be plot
#
# @param choice the graph to plot possible values are "screeplot","quantif","indiv","cor","modalities","mixed","squared loadings". See Details.
# @param level.scale vector(length p) giving the nature of each variable. Possible values: "nom", "ord", "num"
# @param level.scale.supp nature of supplementary variables
# @param supp.var if the user asked to plot the supplementary variable or not
# @param comp axes to use in the plots
# @param nb.comp number of components in the PCAOS model
#
# @return stop function if argument are wrong
#
check.plot.arg <- function(choice,level.scale,level.scale.supp,supp.var, comp, nb.comp,rank){
  if (!(choice %in% c("screeplot","quantif","ind","numeric","qualitative","mixed"))){
    stop("Values of choice should be one of : screeplot,numeric,qualitative,quantif,ind,mixed")
  }
  if (choice == "qualitative"){
    if (!any(level.scale == "nom"|level.scale == "ord")){
      stop("No variable defined as qualitative (i.e nominal or ordinal)")
    }
  }

  if (choice == "ind"){
    if(supp.var == TRUE & !any(level.scale.supp == "nom" | level.scale.supp == "ord")){
      stop("Supplementary variable is not qualitative")
    }
  }

  if (choice == "quantif" & rank == "no.restriction"){
    stop("Quantification plot not available for multiple quantification")
  }

  if (choice == "numeric"){
    if (!any(level.scale == "num")){
      stop("No variable defined as numeric")
    }
    if(supp.var == TRUE & !any(level.scale.supp == "num")){
      stop("Supplementary variable is not numeric")
    }
  }

  # if (choice == "mixed"){
  #   if ((!any(nature == "nom") | !(any(nature == "ord"))) & !any(nature == "num") ){
  #     stop("Variables are not mixed")
  #   }
  # }

  if (any(comp > nb.comp)){
    stop("Wrong axes to plot")
  }

}
