# ordinal.quantification
#
# ordinal.quantification
#
# @param var raw ordinal variables to be quantified
#
# @param t components scores obtain trough PCAOS method
#
# @return
#  \itemize{
#   \item var.quant : optimally quantified variables
#   \item w : loadings
#   \item Yj : quantification of the categories
#   \item Yjhat : rank one quantification of the categories
# }
#
# @details
# Do NOT use this function unless you are ME, a package developer, or a jedi user who really knows what is doing.
#
ordinal.quantification <- function (var, t) {
  decreased <- FALSE
  #missing data
  na <- which(is.na(var))

  nbindiv <- nrow(var)
  var.dis <- dataDIS1 <- dataDIS2 <- dummy.coding(var)$data
  f12<-f121<-f122<-sqrt(diag(t(var.dis)%*%var.dis))
  Yj<-Yj1<-Yj2 <- solve(t(var.dis) %*% var.dis) %*% t(var.dis) %*% t
  ressvd=svd(diag(f12)%*%Yj,nu=1,nv=1)
  qj<-qj1<-qj2<-diag(1/f12)%*%ressvd$u
  rownames(qj) <- rownames(qj1) <- rownames(qj2) <- colnames(dataDIS2)
  aj<-aj1<-aj2<-ressvd$v*ressvd$d[1]
  var.quant.nom <-var.dis %*% qj
  coefnorm<-sqrt(nbindiv)

  #Si quantif bien ordonnee on renvoie les resultats de la quantif nominal
  if(sum(order(qj)==seq(1,ncol(var.dis)))==ncol(var.dis)){
    Yjhat=qj%*%t(aj)
    var.quant <- var.quant.nom* coefnorm
    w<-aj/coefnorm
    return(list(var.quant=var.quant, w=t(w),Yj=Yj,Yjhat=Yjhat))
  }else{
    #Sinon PAVA
    #1.Decreasing
    nbcolDISini1 <- ncol(dataDIS1)
    i=1
    while (i<ncol(dataDIS1)) {
      if (qj1[i] > qj1[i + 1]){
        dataDIS1[, i + 1] <- dataDIS1[,i] + dataDIS1[, i + 1]
        dataDIS1 <- dataDIS1[,-i,drop = FALSE]
        f121<-sqrt(diag(t(dataDIS1)%*%dataDIS1))
        Yj1 <- solve(t(dataDIS1) %*% dataDIS1) %*% t(dataDIS1) %*% t
        if (length(f121)==1) {
          qj1=f121
          aj1=Yj1
        } else {
          ressvd=svd(diag(f121, nrow = length(f121))%*%Yj1,nu=1,nv=1)
          qj1<-diag(1/f121, nrow = length(f121))%*%ressvd$u
          rownames(qj1) <- colnames(dataDIS1)
          aj1<-ressvd$v*ressvd$d[1]
          if (ncol(t)==1) {
            s=as.numeric(sign(aj1*aj))
          } else {
            s=as.numeric(sign(aj1[1]*aj[1]))
          }
          aj1=aj1*s
          qj1=qj1*s
        }
      } else {
        i=i+1
      }
    }
    quantifRESTRICTED1 <- dataDIS1 %*% qj1
    quantifRESTRICTED1[na] <- mean(quantifRESTRICTED1,na.rm = TRUE) #*** ?

    #2.Increasing
    nbcolDISini2 <- ncol(dataDIS2)
    i=1
    while (i<ncol(dataDIS2)) {
      if (qj2[i] < qj2[i + 1]){
        dataDIS2[, i + 1] <- dataDIS2[,i] + dataDIS2[, i + 1]
        dataDIS2 <- dataDIS2[,-i,drop = FALSE]
        f122<-sqrt(diag(t(dataDIS2)%*%dataDIS2))
        Yj2 <- solve(t(dataDIS2) %*% dataDIS2) %*% t(dataDIS2) %*% t
        if (length(f122)==1) {
          qj2=f122
          aj2=Yj2
        } else {
          ressvd=svd(diag(f122, nrow = length(f122))%*%Yj2,nu=1,nv=1)
          qj2<-diag(1/f122, nrow = length(f122))%*%ressvd$u
          aj2<-ressvd$v*ressvd$d[1]
          if (ncol(t)==1) {
            s=as.numeric(sign(aj2*aj))
          } else {
            s=as.numeric(sign(aj2[1]*aj[1]))
          }
          aj2=aj2*s
          qj2=qj2*s
        }
      } else {
        i=i+1
      }
    }

    quantifRESTRICTED2 <- dataDIS2 %*% qj2
    quantifRESTRICTED2[na] <- mean(quantifRESTRICTED2,na.rm = TRUE)  #*** ?

    SSQ1=as.numeric((nbindiv-1)*var(quantifRESTRICTED1))
    SSQ2=as.numeric((nbindiv-1)*var(quantifRESTRICTED2))

    if (abs(SSQ1-SSQ2)<1e-2) {
      if (cor(var.quant.nom,quantifRESTRICTED1)>cor(var.quant.nom,quantifRESTRICTED2)) {
        dataDIS <- dataDIS1
        Yj=Yj1
        aj<-aj1
        qj<-qj1
        var.quant <- (quantifRESTRICTED1)
      } else {
        dataDIS <- dataDIS2
        Yj=Yj2
        aj<-aj2
        qj<-qj2
        var.quant <- (quantifRESTRICTED2)
        decreased <- TRUE
      }
    } else {
      if (SSQ1>SSQ2) {
        dataDIS <- dataDIS1
        Yj=Yj1
        aj<-aj1
        qj<-qj1
        var.quant <- (quantifRESTRICTED1)
      } else {
        dataDIS <- dataDIS2
        Yj=Yj2
        aj<-aj2
        qj<-qj2
        var.quant <- (quantifRESTRICTED2)
      }
    }

    unique.var <- unique(var[,1])
    if (any(is.na(var))){unique.var <- unique(var)[-which(is.na(unique(var))),]}

    if (ncol(dataDIS)<ncol(var.dis)) {
      qjreduced<-qj
      cpt=0
      for (val in sort(unique.var)) {
        cpt=cpt+1
        qj[cpt]<-var.quant[which(var[,1]==val)][1]
      }
      corresp<-rep(0,length(unique.var))
      cpt=0
      for (val in qj) {
        cpt=cpt+1
        corresp[cpt]<-which(qjreduced==val)
      }
      Yjreduced=Yj
      Yj=NULL
      for (d in 1:ncol(Yjreduced)) Yj<-cbind(Yj,Yjreduced[,d][corresp])
    }

    Yjhat= qj%*%t(aj)
    var.quant <- var.quant* coefnorm
    w<-aj/coefnorm

  }

  # #Correction avec gpava si la quantification n'est tjrs pas parfaite
  # if(!sum(order(qj)==seq(1,ncol(var.dis)))==ncol(var.dis)){
  #   if(decreased == TRUE){
  #     qj.order <- gpava(z = length(qj):1, y = qj, weights = 1/f12)$x
  #   }
  #   if(decreased == FALSE){
  #     qj.order <- gpava(z = 1:length(qj), y = qj, weights = 1/f12)$x
  #   }
  #   qqj <- qj.order/as.numeric(sqrt(t(qj.order)%*%qj.order))
  #   #qqj=diag(1/f12)%*%qj.order   #t(qqj)%*%t(var.dis)%*%var.dis%*%qqj=1
  #   var.quant <- var.dis %*% qqj
  #   # centré var.quant
  #   if(qj.order[1] != qj.order[length(qj.order)]){
  #     var.quant=scale(var.quant,center=TRUE,scale=FALSE)
  #   }
  #   var.quant=var.quant/as.numeric(sqrt(t(var.quant)%*%var.quant))
  #   aj=aj%*%t(var.quant.nom)%*%var.quant
  #   Yjhat=qj%*%t(aj)
  #   coefnorm<-sqrt(nbindiv)
  #   var.quant=var.quant*coefnorm
  #   w<-aj/coefnorm
  # }


  return(list(var.quant=var.quant, w=t(w),Yj=Yj,Yjhat=Yjhat))

}

