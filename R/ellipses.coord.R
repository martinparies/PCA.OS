ellipses.coord <- function (var.quali, components, modal,level.conf){
  #products <- unique(var.quali)
  modalites <- levels(as.factor(var.quali))
  nb.modal <- length(modalites)
  coord.ellipse <- list(NULL)

  #Boucle sur chaques modalites de la variable quali
  for (modal in 1:nb.modal){
    #Coordonnees individus pour modal sur CP1
    CP1 <- components[which(var.quali == modalites[modal]),1]
    #Coordonnees individus pour modal sur CP2
    CP2 <- components[which(var.quali == modalites[modal]),2]
    #Moyenne (place le centre de l'ellipse)
    moy = colMeans(components[which(var.quali == modalites[modal]), ])
    #Ecar-type des deux composantes
    sd.CP1 <- sd(CP1)/sqrt(length(CP1))
    sd.CP2 <- sd(CP2)/sqrt(length(CP2))
    #Cov entre CP1 et CP2
    dat <- data.frame(x = CP1, y = CP2)
    dat.cov <- cov(dat)
    dat.cov - dat.cov / nrow(dat)
    #Utilisation fct ellipse du package ellipse
    ellips.coord <- ellipse::ellipse(dat.cov, centre = moy,scale = c(sd.CP1,sd.CP2),level = level.conf,npoints  = 100)
    coord.ellipse[[modal]] <- ellips.coord
  }
  names(coord.ellipse) <- modalites
  return(coord.ellipse)
}
