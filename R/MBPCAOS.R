#' MultiBlock Principal Components Analysis with Optimal Scaling features
#'
#' Perform MBPCAOS
#'
#' @param data a data frame with n rows (individuals) and p columns (numeric, nominal and/or ordinal variables)
#'
#' @param blocks vector(length k) with number of variables in each bloc
#'
#' @param blocks.name vector(length k) with names of each bloc
#'
#' @param level.scale vector(length p) giving the nature of each variable. Possible values: "nom", "ord", "num"
#'
#' @param block.scaling scaling applied to each block. Possible value are : \itemize{
#'   \item "inertia"(default): each quantified block is divided by its total inertia (sum of square).
#'   \item "lambda1" : each quantified block is divided by its the first singular value.
#'   \item "null" : no scaling is applied
#' }
#'
#' @param nb.comp number of components of the model (by default 2)
#'
#' @param maxiter maximum number of iterations.
#'
#' @param threshold the threshold for assessing convergence
#'
#' @param supp.var a vector indicating the indexes of the supplementary variables
#'
#' @param print boolean (TRUE by default), if TRUE ther order of the categories of ordinal variables are print
#'
#' @param init  Intitialization strategy, possible values are :
#' \itemize{
#'   \item "rdm" (default) : random initialisation
#'   \item "svd": components are initialized with the singular value decomposition of the concatenated and pretreated variables (numeric variables are standardized and categorical variables are coded as pseudo disjunctive and weighted by the frequencies of the categories)
#' }

#' @return
#'
#' Dimension reduction
#'  \itemize{
#'   \item weigths : list of weights of the variables (loadings and weights are the same in PCA-like model)
#'   \item components : data.frame with individuals scores for each dimension
#'   \item inertia : percentage and cumulative percentage of variance of the quantified variables explained
#'   }
#' Quantifications
#'   \itemize{
#'   \item quantified.data : optimally quantified variables
#'   \item quantification.categories.nom : list of optimally quantified categories (nominal variables)
#'   \item quantification.categories.ord : list of optimally quantified categories (ordinal variables)
#'   \item level.scale : nature of scaling choosen for each variable
#'   \item data : orginal dataset
#'   }
#'
#'   Blocks results
#' \itemize{
#'   \item block.components : components associated with each block
#'   \item block.weight : weights associated with each block
#'   \item blocks : number of variable in each block
#'   \item blocks.name : name of each block
#' }
#'
#'  Algorithm
#' \itemize{
#'   \item summary : summary of the number of variables according to their nature
#'   \item loss.tot : global loss for all variables
#'   \item stockiter : evolution of the criterion for each ieration
#' }
#'
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
#'data('antibiotic')
#'antb.uses <- antibiotic[,c('Atb.conso','Atb.Sys')]
#'health <- antibiotic[,c('Age','Loss')]
#'vet.practices <- antibiotic[,c(6:15)]
#'antibiotic <- data.frame(antb.uses,health,vet.practices)


#'# Defining blocks
#'blocks.name =  c("antibiotic.uses","Health.of.turkeys","Veterinary.practices")
#'blocks <- c(2,2,10)
#'
#'# Level of scaling
#'level.scale <- rep(NA,ncol(antibiotic))
#'res.nature <- nature.variables(antibiotic)
#'level.scale [res.nature$p.numeric] <- "num"
#'level.scale [res.nature$p.quali] <- "nom"
#'#Warning; the ordinal nature of variables can not be detected automaticaly.
#'level.scale[c(1,14)] <- "ord"
#'
#' # MBPCAOS
#'res.MBPCAOS <- MBPCAOS(data = antibiotic,
#'                      level.scale = level.scale,
#'                       blocks = blocks,
#'                       blocks.name = blocks.name,
#'                       nb.comp = 3)
#'
#'# Blocks graphs
#'plot.MBPCAOS(x = res.MBPCAOS,choice = 'blocks')
#'
#'
#'
#' @export MBPCAOS
#' @export
#'
MBPCAOS <- function(data,
                    level.scale = rep("num",ncol(data)),
                    blocks,
                    blocks.name = paste("Bloc",1:length(blocks)),
                    block.scaling = "inertia",
                    nb.comp = 2,
                    maxiter = 100,
                    threshold = 1e-6,
                    supp.var = NULL,
                    print = TRUE,
                    init = "rdm") {
  #Forcing arguments
  rank.restriction = "one"
  D = 1

  #Checking arguments
  check.arg.MB(data,level.scale,rank.restriction,blocks,blocks.name,print = print)

  #blocks before supp var
  nb.bloc <- length(blocks)
  blocks.list <- list(NULL)
  cumul <- c(1,sapply(1:nb.bloc,function(j){sum(blocks[1:j])}))
  cumul2 <- c(1,sapply(2:length(cumul),function(x)cumul[x]+1))
  blocks.list <- sapply(1:nb.bloc,function(j){cumul2[j]:cumul[j+1]},simplify = F)
  data.list <- sapply(1:nb.bloc,function(j)data[,blocks.list[[j]],drop = F],simplify = F)

  #Supplementary variable
  level.scale.supp = data.var.supp = NULL
  if(!is.null(supp.var)){
    data.var.supp <- data[, supp.var,drop = FALSE]
    level.scale.supp <- level.scale[supp.var]
    data <- data[, -supp.var]
    blocks.supp <- NULL
    compteur =  1
    for (supp in supp.var){
      for (b in 1:nb.block){
        if(supp %in% blocks.list[[b]]){
          blocks.supp[compteur] <- b
          compteur =  compteur + 1
        }
      }
    }

    for(b in blocks.supp){blocks[b] <- blocks[b] - 1}
    level.scale <- level.scale[-supp.var]
    tri.supp <- list()
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

  #blocks after supp var
  nb.block <- length(blocks)
  blocks.list <- list(NULL)
  cumul <- c(1,sapply(1:nb.block,function(j){sum(blocks[1:j])}))
  cumul2 <- c(1,sapply(2:length(cumul),function(x)cumul[x]+1))
  blocks.list <- sapply(1:nb.block,function(j){cumul2[j]:cumul[j+1]},simplify = FALSE)
  data.list <- sapply(1:nb.block,function(j)data[,blocks.list[[j]],drop = F],simplify = FALSE)

  # 1. Nature of active variables
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
  components <- matrix(NA,nb.indiv, nb.comp)
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
  block.components <- list(NULL)
  block.weight <- list(NULL)
  block.explained <- matrix(NA,nb.block,nb.comp)
  b.scale.stock <- list(NULL)

  # 3. Rank one restriction
  rank.restriction.nom <- NULL
  if (rank.restriction == "no.restriction" & any(level.scale == "nom")){
    rank.restriction <- rep("one", nb.var.init)
    rank.restriction[which(level.scale == "nom")] <- "no.restriction"
    rank.restriction.nom <- rep("no.restriction",tri$nbvarNOM)
  } else{
    rank.restriction <- rep("one", nb.var.init)
    if(any(level.scale == "nom")){rank.restriction.nom <- rep("one",tri$nbvarNOM)}
  }

  #4. Initialisation
  # Random
  if(init == "rdm"){
    t<-matrix(rnorm(n=nb.indiv*nb.var.init,mean=0,sd=1),nb.indiv,nb.var.init)
    t<-scale(t)
    ressvd<- svd(t)
  }

  # with svd solution
  if (init == "svd"){
    # #-NUMERIQUE
    dataQUANTI <- NA
    if(!is.null(tri$nbvarNUM)){dataQUANTI <- data.frame(scale(tri$data.num, scale = TRUE,center = TRUE))}
    #-NOMINAL
    dis.nom <- NA
    if(!is.null(tri$nbvarNOM)){
      dis.nom <- dummy.coding(data.frame(data[,tri$emplacement.nom]))$data
      dis.nom <- freq.weight(dis.nom)
    }
    #-ORDINAL
    dis.ord <- NA
    if(!is.null(tri$nbvarORD)){
      dis.ord <- dummy.coding(data.frame(data[,tri$emplacement.ord]))$data
      dis.ord <- freq.weight(dis.ord)
    }

    #Scores prend une solution de SVD
    data.init <- cbind(dataQUANTI,dis.nom,dis.ord)
    col.has.na <- apply(data.init, 2, function(x){any(is.na(x))})
    if (any(col.has.na == TRUE)){data.init <- data.init[,-which(col.has.na)]}
    ressvd<- svd(data.init)
  }

  t<-as.matrix(ressvd$u[,1:nb.comp]*sqrt(nb.indiv))
  t <- scale(t)
  loss = 1000
  loss.by.var = rep(loss/nb.var.init,nb.var.init)
  deltaloss = 10000
  iter.b <- 1

  # 5. OPTIMAL SCALING / MODEL ESTIMATION
  while(continue){

    # 5.1 QUANTIFICATION
    # -Nominal variables
    if (!is.null(tri$nbvarNOM)){
      for (j in 1:tri$nbvarNOM){
        resnomquant<- nominal.quantification(tri$data.nom[,j,drop = FALSE],t,rank.restriction = rank.restriction.nom[j])
        quantified.data[[tri$emplacement.nom[j]]] <- resnomquant$var.quant
        weights[[tri$emplacement.nom[j]]] <- resnomquant$w
        stock.phi[[tri$emplacement.nom[j]]]<- resnomquant$Yj
        stock.phihat[[tri$emplacement.nom[j]]]<-resnomquant$Yjhat
      }
    }
    # -Ordinal variables
    if(!is.null(tri$nbvarORD)){
      for (j in 1:tri$nbvarORD){
        resordquant <- ordinal.quantification(tri$data.ord[,j,drop = FALSE],t)
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
    #Quantified blocks
    #as data.frame
    if  (!any(rank.restriction == "no.restriction")){
      quantified.data.df <- as.data.frame(quantified.data)
      colnames(quantified.data.df)=colnames(data)
      quantified.data.list <- list(NULL)
      #as list
      for (bloc in 1:nb.block){
        quantified.data.list[[bloc]] <- quantified.data.df[,blocks.list[[bloc]],drop = F]
        colnames(quantified.data.list[[bloc]]) <- colnames(data[,blocks.list[[bloc]],drop = FALSE])
      }
    }
    # SCALING COEFF FOR EACH BLOCK
    if (block.scaling == "inertia")   b.scale<- blocks #* nb.indiv
    if (block.scaling == "lambda1")   b.scale <- sapply(1:nb.block,function(j){svd(quantified.data.list[[j]])$d[1]^2})
    if (block.scaling == "null")      b.scale<- rep(nb.indiv,nb.block)                                                      # **********************************  d au carr? *********
    b.scale.stock[[iter.b]] <- b.scale
    iter.b <- iter.b + 1

    # 5. 2 COMPUTATION OF CRITERION OF Optimal Scaling
    lossANCIEN <- loss
    loss.by.var.ancien <- loss.by.var
    deltaloss.ancien<-deltaloss
    if (!any(rank.restriction == 'no.restriction')){
      quantified.data.crit <- quantified.data.df
    }else{
      quantified.data.crit <- quantified.data
    }
    crit = calculcrit.MB(data,level.scale,t,quantified.data.crit,weights,stock.phi, stock.phihat, nb.block,blocks.list,b.scale,D)
    loss.by.var=crit$critd[,1]
    loss=crit$loss/nb.comp
    multipleloss=crit$multipleloss/nb.comp
    singleloss=crit$singleloss/nb.comp
    stockiter=rbind(stockiter, c( iter, loss, multipleloss,singleloss))

    # 5.3 Incrementation and convergence test
    iter <- iter + 1
    # print(iter)
    # print(round(c(loss, multipleloss,singleloss),4))
    deltaloss = (lossANCIEN - loss)
    if (iter>1) stock.delta.loss[iter-1] <- deltaloss
    if ((abs(deltaloss) < threshold) | (iter == maxiter)) {
      continue <- FALSE
      # par(mfrow=c(1, 2))
      # matplot(cbind(loss.by.var.ancien,loss.by.var))
      # plot(stock.delta.loss,main = paste("Delta loss (iteration number",iter-1,")"))
      # print(iter-1)
      # print(round(c( loss, multipleloss,singleloss),4))
      # par(mfrow=c(1, 1))
    }

    #5.4 Computation of components
    #T : eigenvector matrix of X = [..| 1/bk Xktilde| ...]
    tANCIEN<-t
    XX=NULL
    for (bloc in 1:nb.block) XX=cbind(XX,as.matrix(quantified.data.list[[bloc]])/sqrt(b.scale[bloc]))
    if (any(rank.restriction == 'one')){
      ressvd<-svd(XX)
      t<-ressvd$u[,1:nb.comp]*sqrt(nb.indiv)
      t <- scale(t)
      valp=ressvd$d^2
    }

    # if(any(rank.restriction == "no.restriction")){
    #   S<-matrix(0,nb.indiv,nb.comp)
    #   # Computing Matj
    #   S<-matrix(0,nb.indiv,nb.comp)
    #   Matj<- list(NULL)
    #   for (j in 1:nb.var.init){
    #     Matj[[j]]<- quantified.data[[j]] %*% weights[[j]]
    #     S<-S+Matj[[j]]
    #   }
    #   t<-scale(S,center=TRUE,scale=FALSE)
    #   ressvd<-svd(t)
    #   torth<-ressvd$u*sqrt(nb.indiv)
    #   t<-torth
    #   valp=ressvd$d^2
    # }

  }
  ############  end while #####################################

  # 6. Solution
  # Global components
  rownames(components) = rownames(data)
  components <- t
  colnames(components) <- paste("t",1:nb.comp,sep="")
  rownames(components) = rownames(data)

  # 7. Computing inertia
  valinertia=valp[1:nb.comp]/sum(valp)
  cumul=cumsum(valinertia)
  inertia <- data.frame(inertia = round(valinertia,4) * 100,
                        cumulative.percentage.of.inertia = round(cumul,4) * 100)

  #Results for each blocks
  #block.components : tk => block.components
  #block.weight : wk => block.weight
  # percentage of explained inertia for each block k: inertia.k => block.explained
  # In=diag(rep(1,nb.indiv))
  # for (bloc in 1:nb.block) {
  #   tk=matrix(0,nb.indiv,nb.comp)
  #   wk=matrix(0,blocs[bloc],nb.comp)
  #   Xksc=as.matrix(quantified.data.list[[bloc]])/sqrt(b.scale[bloc])
  #   Xx=Xksc
  #   for (h in 1:nb.comp) {
  #     tk[,h]=Xx%*%t(Xx)%*%t[,h]
  #     wk[,h]=t(Xx)%*%t[,h]
  #     Xx = (In-(t[,h]%*%t(t[,h])/nb.indiv))%*%Xx  # residus
  #     #E.k = Xx - (Xx%*%t(Xx)%*%t[,h] %*% t(t[,h,drop = F]) %*% Xx )
  #   }
  #   block.components[[bloc]]=tk
  #   block.weight[[bloc]]=wk
  #   rownames(block.weight[[bloc]]) = colnames(Xksc)
  #   inertia.k=diag(t(t)%*%Xksc%*%t(Xksc)%*%t)/nb.indiv
  #   block.explained[bloc,]=round(inertia.k,4)*100
  # }

  # 8. Results for each blocks
  #block.components : tk => block.components
  #block.weight : wk => block.weight
  # percentage of explained inertia for each block k: inertia.k => block.explained
  for (bloc in 1:nb.block) {
    Xtildetilde = as.matrix(quantified.data.list[[bloc]])/sqrt(b.scale[bloc])
    #tk
    block.components[[bloc]] <- matrix(0,nb.indiv,nb.comp)
    block.components[[bloc]] <- Xtildetilde %*% t(Xtildetilde) %*% t

    #wk
    block.weight[[bloc]] <- matrix(0,blocks[bloc],nb.comp)
    block.weight[[bloc]] <- t(Xtildetilde) %*% t

    #inertia
    inertia.k=diag(t(t)%*%Xtildetilde%*%t(Xtildetilde)%*%t)/nb.indiv
    block.explained[bloc,]=round(inertia.k,4)*100
  }

  names(block.components) = names(block.weight) = blocks.name
  colnames(block.explained) = paste('CP',1:nb.comp,sep = '')
  rownames(block.explained) = blocks.name
  names(weights) <- colnames(data)

  #Contributions of blocks
  #cor.bloc <- data.frame(t(sapply(1:length(res.MBPCAOS$block.components),function(b)diag(cor(res.MBPCAOS$block.components[[b]],res.MBPCAOS$components)))))
  index.group <- unlist(sapply(1:nb.block,function(i) rep(i,blocks[i])))
  contrib <- t(data.frame(weights))
  for (i in 1:ncol(contrib)){
    contrib[,i] <- (contrib[,i]^2) #/ sum(contrib[,i]^2)*100
  }
  contrib.list <- sapply(1:nb.block,function(j) (contrib[blocks.list[[j]],]),simplify = F) #/ (b.scale[j])

  contrib.b <- matrix(NA,nb.block,nb.comp)
  for (i in which(unlist(lapply(contrib.list,is.matrix)))) {
    contrib.list[[i]] <- apply(contrib.list[[i]],2,sum)
  }

  contrib.b <- data.frame(t(data.frame(contrib.list)))
  rownames(contrib.b) <- blocks.name
  colnames(contrib.b) <- paste("CP",1:ncol(contrib.b),sep="")

  # 9. Various results
  stockiter<-stockiter[-1,]
  #colnames(stockiter)=c("iter","loss","multipleloss","singleloss")

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

  # 10. Sign of loadings
  if (any(level.scale == "num")){
    for (j in 1:tri$nbvarNUM) {
      var.quant <- quantified.data[[tri$emplacement.num[j]]]
      var<-data[,tri$emplacement.num[j]]
      w<-weights[[tri$emplacement.num[j]]]
      s=as.numeric(sign(cor(var,var.quant)))
      quantified.data.list[[tri$emplacement.num[j]]] <- var.quant * s
      weights[[tri$emplacement.num[j]]]<-w*s
    }
  }

  if (any(level.scale == "ord")){
    for (j in 1:tri$nbvarORD) {
      var.quant<-quantified.data[[tri$emplacement.ord[j]]]
      #var<-as.numeric(data[,tri$emplacement.ord[j]])
      var <- var.lev <- as.character(data[,tri$emplacement.ord[j]])
      lev <- levels(factor(data[,tri$emplacement.ord[j]]))
      num <- 1
      for (cat in 1:length(lev)){
        var.lev[which(var == lev[cat])] <- num
        num <- num + 1
      }
      w<-weights[[tri$emplacement.ord[j]]]
      s=as.numeric(sign(cor(as.numeric(var.lev),var.quant)))
      quantified.data[[tri$emplacement.ord[j]]] <- var.quant * s
      weights[[tri$emplacement.ord[j]]]<- w*s
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
      rownames(quant.MODAL.nom[[i]]) <- levels(var)
      colnames(quant.MODAL.nom[[i]]) <- paste("CP",1:ncol(var.quant),sep="")
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

  if(!any(rank.restriction == "no.restriction")){
    quantified.data <- matrix(unlist(quantified.data),nrow=nb.indiv,ncol=nb.var.init)
    colnames(quantified.data)= rownames(crit$critd) = colnames(data)
  }

  # 10. Print end message
  if (print == TRUE){
    print(
      paste(
        'MBPCAOS algorithm converged in',
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

  dimension.reduction <- list(weights = lapply(weights,as.vector),
                              components = components,
                              inertia = inertia)

  quantification <- list(quantified.data = quantified.data,
                         quantified.data.list = quantified.data.list,
                         quantification.categories.nom = quant.MODAL.nom,
                         quantification.categories.ord = quant.MODAL.ord,
                         data = data,
                         level.scale = level.scale,
                         summary = summary)

  algo <- list(loss.tot = loss,
               stockiter = stockiter,
               crit.var = crit$critd)

  supp.var <- list(var.supp = data.var.supp,
                   level.scale.supp = level.scale.supp,
                   coord.supp.num = coord.supp.num,
                   coord.supp.quali = coord.supp.quali)

  blocks <- list(block.components = block.components,
                 block.weight = block.weight,
                 contrib.blocks= contrib.b,
                 block.scaling = b.scale,
                 blocks = blocks,
                 blocks.name = blocks.name)

  res <- list(Dimension.reduction = dimension.reduction,
              Quantification = quantification,
              Blocks = blocks,
              Algo = algo,
              Supp.var = supp.var)

  # res <-
  #   list(
  #     weights = lapply(weights,as.vector),
  #     components = components,
  #     quantified.data = quantified.data,
  #     quantified.data.list = quantified.data.list,
  #     block.components = block.components,
  #     block.weight = block.weight,
  #     summary = summary,
  #     quantification.categories.nom = quant.MODAL.nom,
  #     quantification.categories.ord = quant.MODAL.ord,
  #     inertia = inertia,
  #     block.explained=block.explained,
  #     loss.tot = loss,
  #     stockiter = stockiter,
  #     data = data,
  #     blocs = blocs,
  #     blocs.name = blocs.name,
  #     level.scale = level.scale,
  #     level.scale.supp = level.scale.supp,
  #     coord.supp.quali = coord.supp.quali,
  #     coord.supp.num = coord.supp.num,
  #     crit.var = crit$critd,
  #     block.scaling = b.scale,
  #     b.scale.stock = b.scale.stock
  #   )

  class(res) = "MBPCAOS"
  return(res)

}



