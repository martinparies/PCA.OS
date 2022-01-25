#' Principal Component Analysis with Optimal Scaling
#'
#' Perform PCAOS
#'
#' @param data a data frame with n rows (individuals) and p columns (numeric, nominal and/or ordinal variables)
#'
#' @param nature vector(length p) giving the nature of each variable. Possible values: "nom", "ord", "num". The order of categories for an ordinal variable is indicated by it's level.

#' @param nb.comp number of components of the model (by default 2)

#' @param maxiter maximum number of iterations
#'
#' @param threshold threshold for assessing convergence
#'
#' @param D degree of the relation between quantified variables and components (only for numeric variables, default = 1)
#
#' @param supp.var a vector indicating the indexes of the supplementary variables
#'
#' @param rank.restriction restriction of the quantification matrix for nominal variable, possible values are :
#' \itemize{
#'   \item "one" (default) : apply a rank-one restriction (i.e single quantification of each variable)
#'   \item "no.restriction": apply no restriction (i.e multiple quantification of each variable)
#' }
#'
#' @param print.order  boolean (TRUE by default), if TRUE ther order of the categories of ordinal variables are print

#' @return
#'
#' \itemize{
#'   \item weigths : list of weights of the variables (loadings and weights are the same in PCA-like model)
#'   \item components : data.frame with individuals scores for each dimension
#'   \item quantified.data : Optimally quantified variables
#'   \item summary : summary of the number of variables according to their nature
#'   \item quantification.categories.nom : list of optimally quantified categories (nominal variables)
#'   \item quantification.categories.ord : list of optimally quantified categories (ordinal variables)
#'   \item inertia : percentage and cumulative percentage of variance of the quantified variables explained
#'   \item loss.tot : global loss for all variables
#'   \item stockiter : evolution of the criterion for each ieration
#'   \item data : orginal dataset
#'   \item nature : nature of scaling choosen for each variable
#'   \item quali.var.supp : supplementary variable
#'
#' }
#'
#' @examples
#' data (antibiotic)
#' nature <- rep(NA,ncol(antibiotic)) #Setting nature argument
#' nature[c(2,3,4)] <- "num"
#' nature[c(1,5,6,7,8,9,10,11,12,13,14,15)] <- "nom"
#' nature[c(1,15)] <- "ord"
#' nature
#'
#'res.PCAOS <- PCA.OS::PCAOS(
#'  data = antibiotic,
#'  nature = nature,
#'  supp.var = 1
#')
#'
#'PCA.OS::plot.PCAOS(
#'  res.PCAOS = res.PCAOS,
#'  choice = "ind",
#'  coloring = antibiotic$Atb.conso)
#'
#' @author
#' \itemize{
#'   \item Martin PARIES (Maintainer: \email{martin.paries@oniris-nantes.fr})
#'   \item Evelyne Vigneau
#'   \item Stephanie Bougeard
#' }
#'
#' @references
#' Paries, Bougeard, Vigneau (2022), Multivariate analysis of Just-About-Right data with optimal scaling approach. Food Quality and Preference (submit)
#'
#' @export
PCAOS <- function(data,
                  nature = rep("num",ncol(data)),
                  nb.comp = 2,
                  maxiter = 100,
                  threshold = 1e-6,
                  D = 1,
                  rank.restriction = "one",
                  supp.var = NULL,
                  print.order = TRUE) {
  #Checking arguments
  check.arg(data,nature,rank.restriction,print.order)

  #Supplementary variable
  nature.supp = NULL
  tri.supp <- list()
  if(!is.null(supp.var)){
    data.var.supp <- data[, supp.var,drop = FALSE]
    nature.supp <- nature[supp.var]
    data <- data[, -supp.var]
    nature <- nature[-supp.var]
    if(any(nature.supp == "ord"|nature.supp == "nom")){
      tri.supp$nb.var.supp.quali <-  length(which(nature.supp == "ord"|nature.supp == "nom"))
      tri.supp$emplacement.var.supp.quali  <- which(nature.supp =="ord"|nature.supp == "nom")
      tri.supp$data.supp.quali <- data.var.supp[,tri.supp$emplacement.var.supp.quali,drop = FALSE]
    }
    if(any(nature.supp == "num")){
      tri.supp$nb.var.supp.num <-  length(which(nature.supp == "num"))
      tri.supp$emplacement.var.supp.num <- which(nature.supp =="num")
      tri.supp$data.supp.num <- data.var.supp[,tri.supp$emplacement.var.supp.num,drop = FALSE]
    }
  }

  # 1. Nature of variables
  tri <- list()
  if(any(nature == "ord")){
    tri$nbvarORD <-  length(which(nature == "ord"))
    tri$emplacement.ord <- which(nature =="ord")
    tri$data.ord <- data[,tri$emplacement.ord,drop = FALSE]
  }
  if(any(nature == "nom")){
    tri$nbvarNOM <-  length(which(nature == "nom"))
    tri$emplacement.nom <- which(nature =="nom")
    tri$data.nom <- data[,tri$emplacement.nom,drop = FALSE]
  }
  if(any(nature == "num")){
    tri$nbvarNUM <-  length(which(nature == "num"))
    tri$emplacement.num <- which(nature =="num")
    tri$data.num <- data[,tri$emplacement.num,drop = FALSE]
  }

  # 2. Table of results
  nb.indiv <- nrow(data)
  nb.var.init <- ncol(data)
  components <- matrix(NA,nb.indiv, nb.comp )
  weights <- list(NULL)
  quantified.data <- list(NULL)
  stock.phi<- list(NULL)
  stock.phihat<- list(NULL)
  valeurpropres <- rep(NA, nb.comp)
  quant.MODAL.nom <- list(NULL)
  quant.MODAL.ord <- list(NULL)
  quant.MODAL.num <- list(NULL)
  dataINITIALISATION <- matrix(NA,nb.indiv,1)
  compteur <- 1
  iter <- 0
  stockiter<-NULL
  continue <- TRUE
  dis.nomP <- NULL
  dis.ordP <- NULL
  stock.delta.loss <- vector()


  # 3. Rank one restriction
  rank.restriction.nom <- NULL
  if (rank.restriction == "no.restriction" & any(nature == "nom")){
    rank.restriction <- rep("one", nb.var.init)
    rank.restriction[which(nature == "nom")] <- "no.restriction"
    rank.restriction.nom <- rep("no.restriction",tri$nbvarNOM)
  } else{
    rank.restriction <- rep("one", nb.var.init)
    if(any(nature == "nom")){rank.restriction.nom <- rep("one",tri$nbvarNOM)}
  }

  # # Not useful if random initialisation of t
  # #-NUMERIQUE
  # if(!is.null(tri$nbvarNUM)){
  #   dataQUANTI <- scale(tri$data.num, scale = TRUE)
  # }  else {
  #   dataQUANTI <- NULL
  # }
  # #-NOMINAL
  # if(!is.null(tri$nbvarNOM)){
  #   dis.nom <- dummy.coding(data.frame(data[,tri$emplacement.nom]))
  #   dis.nomP <- freq.weight(dis.nom$data)
  # }
  # #-ORDINAL
  # if(!is.null(tri$nbvarORD)){
  #   dis.ord <- dummy.coding(data.frame(data[,tri$emplacement.ord]))
  #   dis.ordP <- freq.weight(dis.ord$data)
  # }

  # 5. Initialisation
  # random
  t<-matrix(rnorm(n=nb.indiv*nb.var.init,mean=0,sd=1),nb.indiv,nb.var.init)
  t<-scale(t)
  ressvd<- svd(t)
  #ressvd<- svd(cbind(dataQUANTI,dis.nomP,dis.ordP))
  t<-as.matrix(ressvd$u[,1:nb.comp]*sqrt(nb.indiv))
  loss=10000
  loss.by.var = rep(loss/nb.var.init,nb.var.init)
  deltaloss=10000

  # 6. OPTIMAL SCALING / MODEL ESTIMATION
  while(continue){
    # 6.1 QUANTIFICATION
    # -Nominal variables
    if (!is.null(tri$nbvarNOM)){
      for (j in 1:tri$nbvarNOM){
        resnomquant<-nominal.quantification(tri$data.nom[,j,drop = FALSE],t,rank.restriction = rank.restriction.nom[j])
        quantified.data[[tri$emplacement.nom[j]]] <- resnomquant$var.quant
        weights[[tri$emplacement.nom[j]]] <- resnomquant$w
        stock.phi[[tri$emplacement.nom[j]]]<- resnomquant$Yj
        stock.phihat[[tri$emplacement.nom[j]]]<-resnomquant$Yjhat
      }
    }
    # -Ordinal variables
    if(!is.null(tri$nbvarORD)){
      for (j in 1:tri$nbvarORD){
        resordquant<-ordinal.quantification(tri$data.ord[,j,drop = FALSE],t)
        quantified.data[[tri$emplacement.ord[j]]]  <- resordquant$var.quant
        weights[[tri$emplacement.ord[j]]] <- resordquant$w
        stock.phi[[tri$emplacement.ord[j]]]<-resordquant$Yj
        stock.phihat[[tri$emplacement.ord[j]]]<-resordquant$Yjhat
      }
    }
    #-Numeric variables
    if (!is.null(tri$nbvarNUM)){
      for (j in 1:tri$nbvarNUM){
        resnumquant<-numeric.quantification(tri$data.num[,j,drop = FALSE],t,D=D)
        quantified.data[[tri$emplacement.num[j]]]  <- resnumquant$var.quant
        weights[[tri$emplacement.num[j]]] <- resnumquant$w
        stock.phi[[tri$emplacement.num[j]]]<-resnumquant$Yj
        stock.phihat[[tri$emplacement.num[j]]]<-resnumquant$Yjhat
      }
    }

    # 6.2  loss function
    lossANCIEN <- loss
    loss.by.var.ancien <- loss.by.var
    deltaloss.ancien<-deltaloss
    critd=matrix(NA,nb.var.init,5)
    colnames(critd)=c("lossv1","lossv2","mloss","sloss","nujh2")
    for (j in 1:nb.var.init) {
      if (nature[j]!="num") Zj=dummy.coding(data[,j,drop=FALSE])$data
      if (nature[j]=="num") {
        var=data[,j,drop=FALSE]
        var=scale(var)*sqrt(nb.indiv/(nb.indiv-1))
        Zj=polynom.data(var,D)
      }
      Yjhat=stock.phihat[[j]]
      Yj=stock.phi[[j]]
      matj<- quantified.data[[j]] %*% weights[[j]]
      critd[j,1]=sum(diag(t(matj - t) %*% (matj - t)))/nb.var.init
      critd[j,2]=sum(diag(t((Zj%*%Yjhat)-t)%*%((Zj%*%Yjhat)-t)))/nb.var.init
      critd[j,3]=sum(diag(t((Zj%*%Yj)-t)%*%((Zj%*%Yj)-t)))/nb.var.init
      critd[j,4]=sum(diag(t(Yjhat-Yj)%*%t(Zj)%*%Zj%*%(Yjhat-Yj)))/nb.var.init
      critd[j,5]=sum(diag((t(Yj)%*%t(Zj)%*%Zj%*%Yj)))/nb.var.init
    }

    critd=critd/(nb.indiv*nb.comp)
    loss=sum(critd[,2])
    multipleloss=sum(critd[,3])
    singleloss=sum(critd[,4])

    stockiter=rbind(stockiter, c( iter, loss, multipleloss,singleloss))
    # print(c( iter, loss, multipleloss,singleloss))
    # print("======================================")
    loss.by.var<-critd[,2]

    # 6.3 Incrementation and convergence test
    iter <- iter + 1
    if (iter > 0) {
      deltaloss = (lossANCIEN - loss)
      if (iter>1) stock.delta.loss[iter] <- deltaloss
      if ((deltaloss < threshold) | (iter == maxiter)) {
        continue <- FALSE
        par(mfrow=c(1, 2))
        # matplot(cbind(loss.by.var.ancien,loss.by.var))
        # plot(stock.delta.loss,main = paste("Delta loss (iteration number",iter,")"))
        # print(iter-1)
        # print(round(c( loss, multipleloss,singleloss),4))
      }
    }

    # 6.4 Computing scores
    tANCIEN<-t
    S<-matrix(0,nb.indiv,nb.comp)
    Matj<- list(NULL)
    for (j in 1:nb.var.init){
      Matj[[j]]<-quantified.data[[j]]%*% weights[[j]]
      S<-S+Matj[[j]]
    }
    t<-scale(S,center=TRUE,scale=FALSE)

    #orthogonalisation of t
    ressvd<-svd(t)
    torth<-ressvd$u*sqrt(nb.indiv)
    t<-torth
  }
  ############  end while #####################################

  # 7. Solution
  names(weights) <- colnames(data)
  rownames(components) = rownames(data)
  components <- t
  colnames(components) <- paste("t",1:nb.comp,sep="")
  rownames(components) = rownames(data)
  colnames(stockiter)=c("iter","loss","multipleloss","singleloss")

  #Supplementary variable
  coord.supp.quali = list(NULL)
  coord.supp.num = list(NULL)
  if(!is.null(supp.var)){
    if(any(nature.supp == "nom" | nature.supp == "ord")){
      for(var in 1:tri.supp$nb.var.supp.quali){
        var.supp.quali <- (tri.supp$data.supp.quali[,var])
        modal <- levels(as.factor(var.supp.quali))
        nb.modal <- length(modal)
        coord.supp.quali[[var]] <- matrix(NA, nb.modal, ncol(components))
        colnames( coord.supp.quali[[var]]) <- colnames(components)
        rownames( coord.supp.quali[[var]]) <- modal
        nbvar <- ncol(data)
        for (mod in 1:nb.modal) {
          coord.supp.quali[[var]][mod,] <- colMeans(components[which(var.supp.quali == modal[mod]), ])
        }
      }
    }

    if(any(nature.supp == "num")){
      for(var in 1:tri.supp$nb.var.supp.num){
        coord.supp.num[[var]] <- cor(scale(tri.supp$data.supp.num),components)
      }
    }
  }

  #Summary
  if(is.null(tri$nbvarNOM)){tri$nbvarNOM <- 0 }
  if(is.null(tri$nbvarORD)){tri$nbvarORD <- 0 }
  if(is.null(tri$nbvarNUM)){tri$nbvarNUM <- 0 }
  summary <- data.frame(NB.var.nom = tri$nbvarNOM,NB.var.ord = tri$nbvarORD,NB.var.num = tri$nbvarNUM)

  if (any(nature == "num")){
    for (j in 1:tri$nbvarNUM) {
      var.quant<-quantified.data[[tri$emplacement.num[j]]]
      var<-data[,tri$emplacement.num[j]]
      w<-weights[[tri$emplacement.num[j]]]
      s=as.numeric(sign(cor(var,var.quant)))
      quantified.data[[tri$emplacement.num[j]]] <- var.quant * s
      weights[[tri$emplacement.num[j]]]<-w*s
    }
  }

  if (any(nature == "ord")){
    for (j in 1:tri$nbvarORD) {
      var.quant<-quantified.data[[tri$emplacement.ord[j]]]
      var<-as.numeric(data[,tri$emplacement.ord[j]])
      w<-weights[[tri$emplacement.ord[j]]]
      s=as.numeric(sign(cor(var,var.quant)))
      quantified.data[[tri$emplacement.ord[j]]] <- var.quant * s
      weights[[tri$emplacement.ord[j]]]<-w*s
    }
  }

  if (any(nature == "nom")){
    for (i in 1:tri$nbvarNOM) {
      var.quant=as.matrix(quantified.data[[tri$emplacement.nom[i]]])
      var=as.factor(data[,tri$emplacement.nom[i]])
      quant.MODAL.nom[[i]]<-matrix(0,nlevels(var),ncol(var.quant))
      for (comp in 1:ncol(var.quant)) {
        for (k in 1:nlevels(var)) {
          quant.MODAL.nom[[i]][k,comp]<-var.quant[which(var==levels(var)[k]),comp][1]
        }
      }
      rownames(quant.MODAL.nom[[i]])<-levels(var)
      colnames(quant.MODAL.nom[[i]])<-paste("CP",1:ncol(var.quant),sep="")
    }
    names(quant.MODAL.nom) <- colnames(data[,nature == "nom"])
  }

  if (any(nature == "ord")){
    for (i in 1:tri$nbvarORD) {
      var.quant=as.matrix(quantified.data[[tri$emplacement.ord[i]]])
      var=as.factor(data[,tri$emplacement.ord[i]])
      quant.MODAL.ord[[i]]<-matrix(0,nlevels(var),ncol(var.quant))
      for (comp in 1:ncol(var.quant)) {
        for (k in 1:nlevels(var)) {
          quant.MODAL.ord[[i]][k,comp]<-var.quant[which(var==levels(var)[k]),comp][1]
        }
      }
      rownames(quant.MODAL.ord[[i]])<-levels(var)
      colnames(quant.MODAL.ord[[i]])<-paste("CP",1:ncol(var.quant),sep="")
    }
    names(quant.MODAL.ord) <- colnames(data[,nature == "ord"])
  }

  if(any(rank.restriction == "one")){
    quantified.data <- matrix(unlist(quantified.data),nrow=nb.indiv,ncol=nb.var.init)
    row.names(quantified.data) <- row.names(data)
    colnames(quantified.data) <- colnames(data)
  }

  # Percentage of inertia for each component
  valinertia=rep(0,nb.comp)
  for (j in 1:nb.var.init) {
    matj=Matj[[j]]
    for (d in 1:nb.comp) {
      valinertia[d]=valinertia[d] + sum(diag(t(matj[,d])%*% matj[,d]))
    }
  }
  valinertia=valinertia/(nb.indiv*nb.var.init)
  cumul=cumsum(valinertia)
  inertia <-
    data.frame(inertia = round(valinertia,4) * 100,
               cumulative.percentage.of.inertia = round(cumul,4) * 100)

  res <-
    list(
      weights = lapply(weights,as.vector),
      components = components,
      quantified.data = quantified.data,
      summary = summary,
      quantification.categories.nom = quant.MODAL.nom,
      quantification.categories.ord = quant.MODAL.ord,
      inertia = inertia,
      loss.tot = loss,
      stockiter = stockiter,
      data = data,
      nature = nature,
      nature.supp = nature.supp,
      coord.supp.quali = coord.supp.quali,
      coord.supp.num = coord.supp.num,
      quali.var.supp = tri.supp$data.supp.quali
    )

  class(res) = "PCAOS"
  return(res)

}
