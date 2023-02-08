# check.plot.MB.arg
#
# check if data and arguments are good to be plot for (MBPCAOS)
#
# @param choice the graph to plot possible values are "screeplot","quantif","indiv","cor","modalities","mixed","squared loadings". See Details.
# @param nature vector(length p) giving the nature of each variable. Possible values: "nom", "ord", "num"
# @param nature.supp nature of supplementary variables
# @param supp.var if the user asked to plot the supplementary variable or not
# @param comp axes to use in the plots
# @param nb.comp number of components in the PCAOS model
#
# @return stop function if argument are wrong
#
check.plot.MB.arg <- function(choice,level.scale,level.scale.supp,supp.var, comp, nb.comp){
  if (!(choice %in% c("screeplot","numeric","qualitative","quantif","ind","all.var","squared.loading","blocs","block.components","block.var"))){
    stop("Values of choice should be one of : screeplot,numeric,qualitative,quantif,ind,mixed,squared.loading, blocs,block.components,block.components")
  }
  if (choice == "qualitative"){
    if (!any(level.scale == "nom")){
      stop("No variable defined as qualitative (i.e nominal or ordinal")
    }
    if(supp.var == TRUE & !any(level.scale.supp == "nom" | level.scale.supp == "ord")){
      stop("Supplementary variable is not qualitative")
    }
  }

  if (choice == "numeric"){
    if (!any(level.scale == "num")){
      stop("No variable defined as numeric")
    }
    if(supp.var == TRUE & !any(level.scale.supp == "num")){
      stop("Supplementary variable is not numeric")
    }
  }

  if (choice == "mixed"){
    if ((!any(level.scale == "nom") | !(any(level.scale == "ord"))) & !any(level.scale == "num") ){
      stop("Variables are not mixed")
    }
  }

  if (any(comp > nb.comp)){
    stop("Wrong axes to plot")
  }

}


