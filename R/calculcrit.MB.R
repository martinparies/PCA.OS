calculcrit.MB <-function (data,level.scale,t,quantified.data,weights,stock.phi, stock.phihat,
                          nb.bloc=NULL,blocs.list=NULL,b.scale=NULL,D) {

  nb.var.init=ncol(data)
  nb.indiv=nrow(data)

  if (is.null(nb.bloc)) {
    nb.bloc=1
    blocs.list=1:nb.var.init
  }
  if (is.null(b.scale)) b.scale=1/nb.indiv

  critd=matrix(NA,nb.var.init,5)
  colnames(critd)=c("lossv1","lossv2","mloss","sloss","nujh2")
  for (bloc in 1:nb.bloc){
    for (j in blocs.list[[bloc]]) {
      if (level.scale[j]!="num") Zj=dummy.coding(data[,j,drop=FALSE])$data
      if (level.scale[j]=="num") {
        var=data[,j,drop=FALSE]
        var=scale(var)*sqrt(nb.indiv/(nb.indiv-1))
        Zj=polynom.data(var,D)
      }
      Yjhat=stock.phihat[[j]]
      Yj=stock.phi[[j]]
      if (is.data.frame(quantified.data)){
        matj <- as.matrix(quantified.data [,j,drop = FALSE]) %*% weights[[j]]
      }else{
        matj <- as.matrix(quantified.data [[j]]) %*% weights[[j]]
      }
      critd[j,1]= sum(diag(t(matj - t) %*% (matj - t)))
      critd[j,2]= sum(diag(t((Zj%*%Yjhat)-t)%*%((Zj%*%Yjhat)-t)))
      critd[j,3]= sum(diag(t((Zj%*%Yj)-t)%*%((Zj%*%Yj)-t)))
      critd[j,4]= sum(diag(t(Yjhat-Yj)%*%t(Zj)%*%Zj%*%(Yjhat-Yj)))
      critd[j,5]= sum(diag((t(Yj)%*%t(Zj)%*%Zj%*%Yj)))
    }
  }
  #critd=critd/(nb.indiv)  # car en divisant par b.scale, on divise par nb.indiv aussi
  loss=0
  multipleloss=0
  singleloss=0
  for (bloc in 1:nb.bloc){
    loss=loss+(sum(critd[blocs.list[[bloc]],1]) / as.numeric(b.scale[bloc]))
    multipleloss=multipleloss+ (sum(critd[blocs.list[[bloc]],3]) / as.numeric(b.scale[bloc]))
    singleloss=singleloss+ (sum(critd[blocs.list[[bloc]],4]) / as.numeric(b.scale[bloc]))
  }

  return(list(critd=critd,loss=loss,multipleloss=multipleloss,singleloss=singleloss))
}
