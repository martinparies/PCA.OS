#' Principal Component Analysis with Optimal Scaling
#'
#' Perform PCAOS
#'
#' @param data a data frame with n rows (individuals) and p columns (numeric, nominal and/or ordinal variables)
#'
#' @param level.scale vector(length p) giving the nature of each variable. Possible values: "nom", "ord", "num". The order of categories for an ordinal variable is indicated by it's level.

#' @param nb.comp number of components of the model (by default 2)

#' @param maxiter maximum number of iterations
#'
#' @param threshold threshold for assessing convergence
#'
#' @param D degree of the relation between quantified variables and components (only for numeric variables, default = 1)
#
#' @param supp.var a vector indicating the indexes of the supplementary variables
#'
#'
#' @param print  boolean (TRUE by default), if TRUE convergence information and order of the categories of ordinal variables are printed.
#'
#' @param init  Intitialization strategy, possible values are :
#' \itemize{
#'   \item "rdm" (default) : random initialisation
#'   \item "svd": components are initialized with the singular value decomposition of the concatenated and pretreated variables (numeric variables are standardized and categorical variables are coded as pseudo disjunctive and weighted by the frequencies of the categories)
#' }
#'
#'
#' @return
#'
#' Dimension reduction
#'  \itemize{
#'   \item weigths : list of weights of the variables (loadings and weights are the same in PCA-like model)
#'   \item components : data.frame with individuals scores for each dimension
#'   \item inertia : percentage and cumulative percentage of variance of the quantified variables explained
#'   }
#'
#' Quantifications
#'   \itemize{
#'   \item quantified.data : optimally quantified variables
#'   \item quantification.categories.nom : list of optimally quantified categories (nominal variables)
#'   \item quantification.categories.ord : list of optimally quantified categories (ordinal variables)
#'   \item level.scale : nature of scaling choosen for each variable
#'   \item data : orginal dataset
#'   }
#'  Algorithm
#' \itemize{
#'   \item summary : summary of the number of variables according to their nature
#'   \item loss.tot : global loss for all variables
#'   \item stockiter : evolution of the criterion for each ieration
#' }
#' Supplementary variables
#' \itemize{
#'   \item var.supp : original supplementary variables
#'   \item level.scale.supp : level of scaling of supplementary variables
#'   \item coord.supp.num : coordinates of supplementary numeric variables (correlation with components)
#'   \item coord.supp.quali : coordinates of qualitatve variables (barycenters)
#' }
#'
#'
#' @examples
#'
#' data (antibiotic)
#'# Level of scaling of each variable
#'# Manually
#'level.scale <- rep(NA,ncol(antibiotic)) #Setting level.scale argument
#'level.scale[c(3,4)] <- "num"
#'level.scale[c(6:14)] <- "nom"
#'level.scale[c(1,15)] <- "ord"

#'# Or using nature.variables()
#'level.scale <- rep(NA,ncol(antibiotic))
#'res.nature <- nature.variables(antibiotic)
#'level.scale [res.nature$p.numeric] <- "num"
#'level.scale [res.nature$p.quali] <- "nom"
#'#Warning; the ordinal nature of variables can not be detected automaticaly.
#'level.scale[c(1,15)] <- "ord"
#'
#'# PCAOS
#'res.PCAOS <- PCAOS(
#' data = antibiotic,
#' level.scale = level.scale,
#' nb.comp = 2)
#'
#'# Plot (individuals)
#'plot.PCAOS(
#'  x = res.PCAOS,
#'  choice = "ind",
#'  coloring.indiv = antibiotic$Atb.conso,
#'  size.legend = 12,
#'  size.label = 4
#')
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
                  level.scale = rep("num",ncol(data)),
                  nb.comp = 2,
                  maxiter = 100,
                  threshold = 1e-6,
                  D = 1,
                  #rank.restriction = "one",
                  supp.var = NULL,
                  print = TRUE,
                  init = 'rdm') {
  #Checking arguments
  rank.restriction = "one"
  check.arg(data,level.scale,rank.restriction,print)

  #Supplementary variable
  level.scale.supp = data.var.supp = NULL
  tri.supp <- list()
  if(!is.null(supp.var)){
    data.var.supp <- data[, supp.var,drop = FALSE]
    level.scale.supp <- level.scale[supp.var]
    data <- data[, -supp.var]
    level.scale <- level.scale[-supp.var]
    if(any(level.scale.supp == "ord"|level.scale.supp == "nom")){
      tri.supp$nb.var.supp.quali <-  length(which(level.scale.supp == "ord"|level.scale.supp == "nom"))
      tri.supp$emplacement.var.supp.quali  <- which(level.scale.supp =="ord"|level.scale.supp == "nom")
      tri.supp$data.supp.quali <- data.var.supp[,tri.supp$emplacement.var.supp.quali,drop = FALSE]
    }
    if(any(level.scale.supp == "num")){
      tri.supp$nb.var.supp.num <-  length(which(level.scale.supp == "num"))
      tri.supp$emplacement.var.supp.num <- which(level.scale.supp =="num")
      tri.supp$data.supp.num <- data.var.supp[,tri.supp$emplacement.var.supp.num,drop = FALSE]
    }
  }

  # 1. Nature of variables
  tri <- list()
  if(any(level.scale == "ord")){
    tri$nbvarORD <-  length(which(level.scale == "ord"))
    tri$emplacement.ord <- which(level.scale =="ord")
    tri$data.ord <- data[,tri$emplacement.ord,drop = FALSE]
  }
  if(any(level.scale == "nom")){
    tri$nbvarNOM <-  length(which(level.scale == "nom"))
    tri$emplacement.nom <- which(level.scale =="nom")
    tri$data.nom <- data[,tri$emplacement.nom,drop = FALSE]
  }
  if(any(level.scale == "num")){
    tri$nbvarNUM <-  length(which(level.scale == "num"))
    tri$emplacement.num <- which(level.scale =="num")
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
  if (rank.restriction == "no.restriction" & any(level.scale == "nom")){
    rank.restriction.nom <- rep("one", nb.var.init)
    rank.restriction.nom[which(level.scale == "nom")] <- "no.restriction"
    rank.restriction.nom <- rep("no.restriction",tri$nbvarNOM)
  } else{
    #rank.restriction <- rep("one", nb.var.init)
    if(any(level.scale == "nom")){rank.restriction.nom <- rep("one",tri$nbvarNOM)}
  }

  # 4. Initialisation
  #by svd solution
  if(init == 'svd'){
    if(!is.null(tri$nbvarNUM)){
      dataQUANTI <- scale(tri$data.num, scale = TRUE)
    }  else {
      dataQUANTI <- NULL
    }
    #-NOMINAL
    if(!is.null(tri$nbvarNOM)){
      dis.nom <- dummy.coding(data.frame(data[,tri$emplacement.nom]))
      dis.nomP <- freq.weight(dis.nom$data)
    }
    #-ORDINAL
    if(!is.null(tri$nbvarORD)){
      dis.ord <- dummy.coding(data.frame(data[,tri$emplacement.ord]))
      dis.ordP <- freq.weight(dis.ord$data)
    }
    init.dat <- cbind(dataQUANTI,dis.nomP,dis.ordP)
    for (i in 1:ncol(init.dat)){init.dat[which(is.na(init.dat[,i])),i] <- 0}
    ressvd<- svd(init.dat)
  }
  # random
  if(init == 'rdm'){
    t<-matrix(rnorm(n=nb.indiv*nb.var.init,mean=0,sd=1),nb.indiv,nb.var.init)
    t<-scale(t)
    ressvd<- svd(t)
  }

  t<-as.matrix(ressvd$u[,1:nb.comp]*sqrt(nb.indiv))
  loss=10000
  loss.by.var = rep(loss/nb.var.init,nb.var.init)
  deltaloss=10000

  # 5. OPTIMAL SCALING / MODEL ESTIMATION
  while(continue){
    # 5.1 QUANTIFICATION
    # -Qualitative variables
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

    # 5.2  loss function
    lossANCIEN <- loss
    loss.by.var.ancien <- loss.by.var
    deltaloss.ancien<-deltaloss
    critd=matrix(NA,nb.var.init,5)
    colnames(critd)=c("lossv1","lossv2","mloss","sloss","nujh2")
    for (j in 1:nb.var.init) {
      if (level.scale[j]!="num") Zj=dummy.coding(data[,j,drop=FALSE])$data
      if (level.scale[j]=="num") {
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

    # 5.3 Incrementation and convergence test
    iter <- iter + 1
    if (iter > 0) {
      deltaloss = (lossANCIEN - loss)
      if (iter>1) stock.delta.loss[iter] <- deltaloss
      if ((deltaloss < threshold) | (iter == maxiter)) {
        continue <- FALSE
        # matplot(cbind(loss.by.var.ancien,loss.by.var))
        # plot(stock.delta.loss,main = paste("Delta loss (iteration number",iter,")"))
        # print(iter-1)
        # print(round(c( loss, multipleloss,singleloss),4))
      }
    }

    #5.4 Computation of components
    if(rank.restriction == "one"){
      XX=NULL
      ressvd<-svd(data.frame(quantified.data))
      t<-ressvd$u[,1:nb.comp]*sqrt(nb.indiv) # pour que t't = NI
      t<-scale(t,center=TRUE,scale=FALSE)
    }
    if(rank.restriction == "no.restriction"){
      S<-matrix(0,nb.indiv,nb.comp)
      # Computing Matj
      S<-matrix(0,nb.indiv,nb.comp)
      Matj<- list(NULL)
      for (j in 1:nb.var.init){
        Matj[[j]]<- quantified.data[[j]] %*% weights[[j]]
        S<-S+Matj[[j]]
      }
      t<-scale(S,center=TRUE,scale=FALSE)
      ressvd<-svd(t)
      torth<-ressvd$u*sqrt(nb.indiv)
      t<-torth

    }

  }
  ############  end while #####################################

  # 6. Solution
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
    if(any(level.scale.supp == "nom" | level.scale.supp == "ord")){
      for(var in 1:tri.supp$nb.var.supp.quali){
        var.supp.quali <- (tri.supp$data.supp.quali[,var])
        modal <- levels(as.factor(var.supp.quali))
        nb.modal <- length(modal)
        coord.supp.quali[[var]] <- matrix(NA, nb.modal, ncol(components))
        colnames( coord.supp.quali[[var]]) <- colnames(components)
        rownames( coord.supp.quali[[var]]) <- modal
        nbvar <- ncol(data)
        for (mod in 1:nb.modal) {
          coord.supp.quali[[var]][mod,] <- colMeans(components[which(var.supp.quali == modal[mod]), ,drop = F])
        }
      }
    }
    names(coord.supp.quali) <- colnames(data.var.supp[,which(level.scale.supp == 'nom'|level.scale.supp =='ord')])

    if(any(level.scale.supp == "num")){
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
  summary <- list(summary = summary,rank = rank.restriction)

  #Contribution of variables
  contrib <-  do.call(rbind.data.frame,weights)
  for (i in 1:nb.comp){contrib[,i] <- contrib[,i]^2}
  for (i in 1:nb.comp){contrib[,i] <- contrib[,i] / sum(contrib[,i])*100}
  rownames(contrib) <- names(weights)
  colnames(contrib) <- paste("CP",1:ncol(contrib),sep="")


  # 7. Sign of loadings
  if (any(level.scale == "num")){
    for (j in 1:tri$nbvarNUM) {
      var.quant<-quantified.data[[tri$emplacement.num[j]]]
      var<-data[,tri$emplacement.num[j]]
      na <- which(is.na(var))
      var[na] <- mean(var,na.rm = T)
      w<-weights[[tri$emplacement.num[j]]]
      s=as.numeric(sign(cor(var,var.quant)))
      quantified.data[[tri$emplacement.num[j]]] <- var.quant * s
      weights[[tri$emplacement.num[j]]]<-w*s
    }
  }

  if (any(level.scale == "ord")){
    for (j in 1:tri$nbvarORD) {
      var.quant <- quantified.data[[tri$emplacement.ord[j]]]
      var <- as.numeric(data[,tri$emplacement.ord[j]])
      na <- which(is.na(var))
      var[na] <- mean(var,na.rm = T)
      w<-weights[[tri$emplacement.ord[j]]]
      s=as.numeric(sign(cor(var,var.quant)))
      quantified.data[[tri$emplacement.ord[j]]] <- var.quant * s
      weights[[tri$emplacement.ord[j]]]<-w*s
    }
  }

  if (any(level.scale == "nom")){
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
    names(quant.MODAL.nom) <- colnames(data[,level.scale == "nom"])
  }

  if (any(level.scale == "ord")){
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
    names(quant.MODAL.ord) <- colnames(data[,level.scale == "ord"])
  }

  # 8. Computing inertia
  S<-matrix(0,nb.indiv,nb.comp)
  Matj<- list(NULL)
  for (j in 1:nb.var.init){
    Matj[[j]]<- quantified.data[[j]] %*% weights[[j]]
    S<-S+Matj[[j]]
  }
  #t<-scale(S,center=TRUE,scale=FALSE)

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

  if(rank.restriction == "one"){
    quantified.data <- matrix(unlist(quantified.data),nrow=nb.indiv,ncol=nb.var.init)
    row.names(quantified.data) <- row.names(data)
    colnames(quantified.data) <- colnames(data)
  }
  if(rank.restriction == "no.restriction"){
    names(quantified.data) <-  colnames(data)
  }


  # res <-
  #   list(
  #     weights = lapply(weights,as.vector),
  #     components = components,
  #     quantified.data = quantified.data,
  #     summary = summary,
  #     quantification.categories.nom = quant.MODAL.nom,
  #     quantification.categories.ord = quant.MODAL.ord,
  #     inertia = inertia,
  #     loss.tot = loss,
  #     stockiter = stockiter,
  #     data = data,
  #     level.scale = level.scale,
  #     level.scale.supp = level.scale.supp,
  #     coord.supp.quali = coord.supp.quali,
  #     coord.supp.num = coord.supp.num,
  #     quali.var.supp = tri.supp$data.supp.quali
  #   )

  # 9 . Final results list
  dimension.reduction <-
    list(
      components = components,
      weights = lapply(weights, as.vector),
      inertia = inertia,
      contrib.var = contrib
    )

  quantification <-
    list(
      quantified.data = quantified.data,
      quantification.categories.nom = quant.MODAL.nom,
      quantification.categories.ord = quant.MODAL.ord,
      data = data,
      level.scale = level.scale,
      summary = summary)

  algo <- list(loss.tot = loss,
               stockiter = stockiter)

  supp.var <- list(var.supp = data.var.supp,
                   level.scale.supp = level.scale.supp,
                   coord.supp.num = coord.supp.num,
                   coord.supp.quali = coord.supp.quali)

  # 10. Print end message
  if (print == TRUE){
    print(
      paste(
        'PCAOS algorithm converged in',
        nrow(stockiter),
        'iteration.'))
    print(
      paste('The',
            nb.comp,
            'optimized components explain a total of',
            inertia[nrow(inertia),2],'% of the',
            ncol(data),'quantified variables'
      )
    )
  }

  res <- list(Dimension.reduction = dimension.reduction,
              Quantification = quantification,
              Algo = algo,
              Supp.var = supp.var)

  class(res) = "PCAOS"
  return(res)

}
