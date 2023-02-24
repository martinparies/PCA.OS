# nominal.quantification
#
# nominal.quantification
#
# @param var raw nominal variables to be quantified
#
# @param t components scores obtain trough PCAOS method
#
# @param rank.restriction whether the quantification should be restricted or not

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
nominal.quantification <- function (var,t,rank.restriction){
  nbindiv <- nrow(var)
  var.dis <- dummy.coding(var)$data
  f12 <- sqrt(diag(t(var.dis)%*%var.dis))
  Yj <- solve(t(var.dis) %*% var.dis) %*% t(var.dis) %*% t
  coefnorm<-sqrt(nbindiv)
  if (rank.restriction == "one") param=1
  if (rank.restriction == "no.restriction") param=min(nrow(Yj), ncol(Yj))
  ressvd=svd(diag(f12)%*%Yj,nu=param,nv=param)
  qj=diag(1/f12)%*%ressvd$u
  if ((rank.restriction == "one")| (param==1))  aj=ressvd$v*ressvd$d[1]
  if ((rank.restriction == "no.restriction") & (param>1))    aj=ressvd$v%*%diag(ressvd$d[1:param])
  Yjhat=qj%*%t(aj)
  var.quant <- var.dis %*% qj
  var.quant <- var.quant * coefnorm
  w <- aj/coefnorm

  return(list(var.quant=var.quant,w=t(w),Yj=Yj,Yjhat=Yjhat,qj=qj))
}
