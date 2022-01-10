#' check.arg
#'
#' check if data and arguments are good to be analysed by PCAOS function
#'
#' @param data raw variables (vector or matrix)
#' @param nature vector(length p) giving the nature of each variable. Possible values: "nom", "ord", "num"
#' @param rank.restriction restriction of the quantification matrix for nominal variable
#'
#' @return stop function if argument or data ar wrong
#'
check.arg <- function(data,nature,rank.restriction){
  check = list(NULL)
  #Structure of data
  if(!is.data.frame(data)){
    stop("Argument data should be a data.frame")
    check$data = TRUE
  }else{
    #nature
    if (ncol(data) != length(nature)){
      stop(paste("Error, the length of nature is different from the number of variable data."))
      check$nature = TRUE
    }

    #Données manquantes
    if(any(is.na(data))){
      var.with.na <- colnames(data[,which(sapply(1:ncol(data),function(i){any(is.na(data[,i]))})),drop=F])
      stop(paste("Error",var.with.na,"have missing values. Please consider imputation.",sep= " "),)
      check$missing = TRUE
    }

    #Restriction
    if(!(rank.restriction %in% c("one","no.restriction"))){
      stop(paste("Error, rank.restriction argument should be one of these: one; no.restriction"))
      check$restriction = TRUE
    }

    #Nature of variables
    mix <- mix.data(data)

    if(any(nature == "num")){
      data.num = data[,which(nature == "num")]
      if(mix$nb.numeric != length(which(nature == "num"))){
        var.facto.def.as.num = colnames(data.num[,which(sapply(1:ncol(data.num),function(i){any(is.factor(data.num[,i]))})),drop=F])
        message(paste("Warning, factors variables ",var.facto.def.as.num ," defined as numeric"))
        check$nature.num = TRUE
      }
    }

    if(any(nature == "nom")){
      data.nom = data[,which(nature == "nom"| nature =="ord")]
      if(mix$nb.quali != length(which(nature == "nom" | nature =="ord"))){
        var.nom.def.as.num = colnames(data.nom[,which(sapply(1:ncol(data.nom),function(i){any(is.numeric(data.nom[,i]))})),drop=F])
        stop(paste("Warning, numeric variables",var.nom.def.as.num ,"defined as factor"))
        check$nature.nom = TRUE
      }
      #A variable with only one categories
      table = apply(data.nom,2,table)
      length =lapply(table,length)
      if(any(unlist(length) == 1 )){
        var.one.cat = names(length[which(unlist(length) == 1 )])
        stop(paste("Error variable",var.one.cat,"has only one categorie. Factor analysis might not be suitable for this variable."))
        check$one.cat = TRUE
      }
    }


  }
}
