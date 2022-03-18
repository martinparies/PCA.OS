# check.arg
#
# check if data and arguments are good to be analysed by PCAOS function
#
# @param data raw variables (vector or matrix)
# @param level.scale vector(length p) giving the nature of each variable. Possible values: "nom", "ord", "num"
# @param rank.restriction restriction of the quantification matrix for nominal variable
# @param print.order  boolean (TRUE by default), if TRUE ther order of the categories of ordinal variables are print
#
# @return stop function if argument or data ar wrong
#
check.arg <- function(data,level.scale,rank.restriction,print.order){
  check = list(NULL)
  #Structure of data
  if(!is.data.frame(data)){
    stop("Argument data should be a data.frame")
  }else{
    #nature
    if (ncol(data) != length(level.scale)){
      stop(paste("Error, the length of level.scale is different from the number of variable data."))
    }

    #Données manquantes
    if(any(is.na(data))){
      var.with.na <- colnames(data[,which(sapply(1:ncol(data),function(i){any(is.na(data[,i]))})),drop=F])
      stop(paste("Error",var.with.na,"have missing values. Please consider imputation.",sep= " "),)
    }

    #Restriction
    if(!(rank.restriction %in% c("one","no.restriction"))){
      stop(paste("Error, rank.restriction argument should be one of these: one; no.restriction"))
    }

    #Nature of variables
    mix <- mix.data(data)

    # if(any(nature == "num")){
    #   data.num = data[,which(nature == "num"),drop = F]
    #   if(mix$nb.numeric != length(which(nature == "num"))){
    #     var.facto.def.as.num = colnames(data.num[,which(sapply(1:ncol(data.num),function(i){any(is.factor(data.num[1,i]))})),drop=F])
    #     #var.character.def.as.num = colnames(data.num[,which(sapply(1:ncol(data.num),function(i){any(is.character(data.num[,i]))})),drop=F])
    #     stop(paste("Warning, factors variables ",var.facto.def.as.num ," defined as numeric"))
    #   }
    # }

    if(any(level.scale == "nom") | any(level.scale == "ord")){
      data.nom <-  data[,which(level.scale == "nom" | level.scale =="ord"),drop = F]
      #A variable with only one categories
      table = sapply(1:ncol(data.nom),function(x){table(data.nom[,x])},simplify = F)
      length = lapply(table,length)
      if(any(unlist(length) == 1 )){
        var.one.cat = names(length[which(unlist(length) == 1 )])
        stop(paste("Error variable",var.one.cat,"has only one categorie. Factor analysis might not be suitable for this variable."))
      }
    }

    if(any(level.scale == "ord") & print.order == TRUE){
      data.ord <-  data[,which(level.scale =="ord"),drop = F]
      order.detect <- list(NULL)
      order.detect <- sapply(1:ncol(data.ord),function(var){levels(as.factor(data.ord[,var]))},simplify = F)
      names(order.detect) <- colnames(data.ord)
      for(var.ord in 1:ncol(data.ord)){
        print(paste("Order detect for variable",colnames(data.ord[,var.ord,drop = F]),"is"))
        print(order.detect[[var.ord]])
      }
    }
  }
}
