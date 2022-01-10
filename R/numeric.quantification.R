#' numeric.quantification
#'
#' numeric.quantification
#'
#' @param var raw numeric variables to be quantified
#' @param D degree assumed between variable and component
#' @param t components scores obtain trough PCAOS method
#'
#' @return
#'  \itemize{
#'   \item var.quant : optimally quantified variables
#'   \item w : loadings
#'   \item Yj : quantification of the categories
#'   \item Yjhat : rank one quantification of the categories
#' }
#'
#' @details
#' Do NOT use this function unless you are ME, a package developer, or a jedi user who really knows what is doing.
#'
numeric.quantification <- function (var,t,D){
nbindiv <- nrow(var)
var=scale(var)*sqrt(nbindiv/(nbindiv-1))
var.poly=polynom.data(var,D)
# var.poly <- rep(1,nbindiv)
# for (puiss in 1:D) {var.poly <- cbind(var.poly,var^puiss)}
# var.poly=as.matrix(var.poly)
Yj<- solve(t(var.poly) %*% var.poly) %*%  t(var.poly) %*% t
ressvd=svd(Yj)
#qj=ressvd$u[,1]
qj=ressvd$u[,1]
aj=ressvd$v[,1]*ressvd$d[1]
Yjhat=qj%*%t(aj)
var.quant <- var.poly %*% qj
w<-aj
#print(sum((var.quant%*%as.vector(aj)-t)^2))
#ne pas faire ce changement de signe avant que l'algo n'ait convergÃ©
#var.quant <- var.quant * as.numeric(sign(cor(var,var.quant)))
# changer le signe des weights aussi


return(list(var.quant=var.quant,w=w,Yj=Yj,Yjhat=Yjhat))
}
