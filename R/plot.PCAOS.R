#' plot.PCAOS
#'
#' Visualisation of results from PCAOS method. See details for available plot.
#'
#' @param x an object of class PCAOS
#' @param choice the available graphs are "screeplot","quantif","ind","numeric","qualitative","all.var". See Details.
#' @param comp a length 2 vector with the components to plot.
#' @param supp.var TRUE or FALSE; if TRUE supplementary variables are added in factorial representation (individuals and variables plot).
#' @param size.label size of label in graphs (all plots).
#' @param size.legend  size of label in graphs (all plots).
#' @param sub.var.quantif a vector with variable of interest (quantification plots).
#' @param coloring.indiv a vector of length N to color individuals. If NULL, no coloring is applied (individuals plot).
#' @param ellipse boolean (FALSE by default), if TRUE, draw ellipses around categories of the qualitative variable considered as supplementary (individuals plot).
#' @param level.conf level of confidence ellipses (individuals plot).
#' @param min.contribution  variables with a contribution (i.e loading) lower than this value will not be plotted in the 'all.var' graph (useful for dataset with a lot of variables) (all.var plot).
#' @param label.cat If == 'var+cat', the name of the variable is included in the labels of the categories; if == 'cat', only the name of the categorie is plotted (name of categories should be unique) (qualitative and all.var plot).
#' @param ordinal.as.direction boolean (FALSE by default); if TRUE ordinal variables are represented as vectors, from the first categorie to the last one (qualitative and all.var plot).
#' @param label.size.freq  boolean (FALSE by default); if TRUE size of categories are proportional to their citation frequencies (qualitative and all.var plot).
#' @param ... further arguments passed to or from other methods, such as cex, cex.main, ...

#' @return
#'A ggplot object
#'
#' @details
#'  \itemize{
#'   \item  screeplot: Representation of the percentage of inertia restituates (Y), for each component (X).
#'   \item  quantif: Reprensetation of the quantification of variables trought Optimal Scaling, with original variables (X) and quantified variables (Y).Possibility to select one or more variables of interest with the argument "var.sub".
#'   \item  ind: factorial representation of individuals
#'   }
#'
#' For numeric variables
#'   \itemize{
#'   \item  numeric: factorial representation of numeric variables (also called loading plot) Each numeric variable is represented by it's weight/loadings
#'   }
#' For qualitative (i.e nominal and ordinal) variables
#' \itemize{
#'   \item  qualitative: factorial representation of qualitatives variables trough the representation of it's categories. Coordinates of each category is calculted such as the single quantification of the category multiplied by the loading of the associated variable (rank.restriction = one). Or by averaging, per component, the principal component scores for all individuals in the same categories of a particular variable (rank.restriction = no.restriction).
#' }
#' For All  variables
#'  \itemize{
#'   \item  all.var : factorial representation of all variables (weight for numeric variables, and categories for qualitative variables)
#' }
#'
#' All graph are ggplot object
#'
#' @examples
#' data (antibiotic)
#' level.scale <- rep(NA,ncol(antibiotic)) #Setting level.scale argument
#' level.scale[c(2,3,4)] <- "num"
#' level.scale[c(1,5,6,7,8,9,10,11,12,13,14,15)] <- "nom"
#' level.scale[c(1,15)] <- "ord"
#' level.scale
#'
#'res.PCAOS <- PCA.OS::PCAOS(
#'  data = antibiotic,
#'  level.scale = level.scale,
#'  supp.var = c(1,2)
#')
#'
#'#Individuals graph
#'PCA.OS::plot.PCAOS(x = res.PCAOS,choice = "ind",coloring = antibiotic$Atb.conso)
#'PCA.OS::plot.PCAOS(x = res.PCAOS,choice = "ind",supp.var = TRUE,ellipse = TRUE)
#'
#'#Quantifications
#'PCA.OS::plot.PCAOS(x = res.PCAOS,choice = "quantif",sub.var.quantif =  c(4,8))
#'
#'#Qualitative variables
#'PCA.OS::plot.PCAOS(x = res.PCAOS,choice = "qualitative",supp.var = TRUE)
#'
#'#Numeric variables
#'PCA.OS::plot.PCAOS(x = res.PCAOS,choice = "numeric",supp.var = TRUE)
#'
#'#All variables
#'PCA.OS::plot.PCAOS(x = res.PCAOS,choice = "all.var",supp.var = TRUE)
#'
#' @author
#' \itemize{
#'   \item Martin PARIES (Maintainer: \email{martin.paries@oniris-nantes.fr})
#'   \item Evelyne Vigneau
#'   \item Stephanie Bougeard
#' }
#'
#' @rdname plot.PCAOS
#' @export plot.PCAOS
#' @export
#'
plot.PCAOS <-
  function(x,
           choice = "ind",
           comp = c(1,2),
           coloring.indiv = NULL,
           supp.var = FALSE,
           sub.var.quantif = NULL,
           ellipse = FALSE,
           level.conf = 0.95,
           size.label = 3.5,
           size.legend = 10,
           min.contribution = 0,
           label.size.freq = FALSE,
           ordinal.as.direction = FALSE,
           label.cat = 'var+cat',...) {

    res.PCAOS <- x
    if(!inherits(res.PCAOS,"PCAOS")) stop("Non convenient object")

    nb.comp <- ncol(res.PCAOS$Dimension.reduction$components)
    check.plot.arg(
      choice,
      res.PCAOS$Quantification$level.scale,
      res.PCAOS$Supp.var$level.scale.supp,
      supp.var,
      comp,
      nb.comp,
      rank = res.PCAOS$Quantification$summary$rank
    )

    #INFORMATION ABOUT THE RESULTS
    level.scale <- res.PCAOS$Quantification$level.scale
    nom.comp <- c(paste("PC",comp[1],sep=""),paste("PC",comp[2],sep=""))
    inertie <- res.PCAOS$Dimension.reduction$inertia

    #Dimension reduction
    components <- res.PCAOS$Dimension.reduction$components
    weights <- res.PCAOS$Dimension.reduction$weights

    #QUANTIFICATIONS
    data <- res.PCAOS$Quantification$data
    data.quantified <- res.PCAOS$Quantification$quantified.data
    variables <- colnames(data)
    nb.var <- length(variables)
    quantification <- vector("list", length = nb.var)
    quantification.nom <- res.PCAOS$Quantification$quantification.categories.nom
    quantification.ord <- res.PCAOS$Quantification$quantification.categories.ord
    var.quali <- which(level.scale == "nom" | level.scale =="ord")
    var.nom <- which(level.scale == "nom" )
    var.ord <- which(level.scale == "ord" )

    compteur <- 1
    for (i in which(level.scale == "nom")) {
      quantification[[i]] <- res.PCAOS$Quantification$quantification.categories.nom[[compteur]]
      compteur <- compteur + 1
    }
    compteur <- 1
    for (i in which(level.scale == "ord")) {
      quantification[[i]] <- res.PCAOS$Quantification$quantification.categories.ord[[compteur]]
      compteur <- compteur + 1
    }
    compteur <- 1
    if(res.PCAOS$Quantification$summary$rank == "one"){
      for (i in which(level.scale == "num")){
        quantification[[i]] <- cbind(data[,i],data.quantified[,i])
        compteur <- compteur + 1
      }
    }
    if(res.PCAOS$Quantification$summary$rank == "no.restriction"){
      for (i in which(level.scale == "num")){
        quantification[[i]] <- cbind(data[,i],data.quantified[[i]])
        compteur <- compteur + 1
      }
    }
    names (quantification) <- variables

    #Levels of scaling
    level.scale <- res.PCAOS$Quantification$level.scale
    level.scale.supp <- res.PCAOS$Quantification$level.scale.supp

    #Var.supp
    var.supp <- res.PCAOS$Supp.var$var.supp
    level.scale.supp <- res.PCAOS$Supp.var$level.scale.supp
    coord.supp.num <- res.PCAOS$Supp.var$coord.supp.num
    coord.supp.quali <- res.PCAOS$Supp.var$coord.supp.quali

    #SCREEPLOT
    if(choice == "screeplot"){
      screeplot <-
        ggplot2::ggplot(inertie, ggplot2::aes(y = inertie[,1], x = (1:nrow(inertie)))) +
        ggplot2::scale_x_discrete(labels = "") +
        ggplot2::geom_bar(stat = "identity", width = 0.5, fill = "steelblue") +
        ggplot2::geom_text(
          ggplot2::aes(label = paste(inertie[,1], "%")),
          vjust = -0.8,
          color = "black",
          size = 3.5
        ) +
        ggplot2::theme_minimal() +
        ggplot2::xlab("Component") +
        ggplot2::ylab("Percentage of inertia restituate (%)") +
        ggplot2::ggtitle("Screeplot") +
        ggplot2::theme_classic(base_size = size.legend)

      return (screeplot)
    }

    #LOADPLOT
    if(choice == "all.var" & !any(level.scale == 'nom'| level.scale == 'ord')){choice <- 'numeric'}
    if(choice == "qualitative" & !any(level.scale == 'nom'| level.scale == 'ord')){choice <- 'numeric'}
    if (choice == "numeric" ){
      weight.num <- data.frame(t(data.frame(weights)))[which(level.scale == "num"),]
      graph.var <-
        ggplot2::ggplot(data = data.frame(weight.num),
                        ggplot2::aes (
                          x = as.numeric(weight.num[,comp[1]]),
                          y = as.numeric(weight.num[,comp[2]])
                        )) +
        ggplot2::geom_label(ggplot2::aes(label = row.names(weight.num)), color = "black",size = size.label) +
        ggplot2::ggtitle("Factorial representation of numeric variables") +
        ggplot2::xlab(paste(nom.comp[1], inertie[comp[1],1]," %")) +
        ggplot2::ylab(paste(nom.comp[2], inertie[comp[2],1]," %")) +
        ggplot2::geom_hline(
          yintercept = 0,
          linetype = "dotted",
          color = "black",
          size = 1
        ) +
        ggplot2::geom_vline(
          xintercept = 0,
          linetype = "dotted",
          color = "black",
          size = 1
        ) +
        ggplot2::theme_classic(base_size = size.legend) + ggplot2::annotate(
          geom = "segment",
          x = rep(0,nrow(weight.num)),
          xend = weight.num[,comp[1]],
          y = rep(0,nrow(weight.num)),
          yend = weight.num[,comp[2]],
          col = "black",
          arrow = ggplot2::arrow(length = grid::unit(0.1, "cm")),size = 0.70
        )  +
        ggplot2::annotate("path",
                      x = 0 + 1 * cos(seq(0, 2 * pi, length.out = 100)),
                      y = 0 + 1* sin(seq(0, 2 * pi, length.out = 100)),size = 0.50)


      #Si il y a une variable quantitative supplementaire
      if (supp.var == TRUE){
        loading.supp <- do.call(rbind.data.frame,coord.supp.num)
        graph.var <-
          graph.var  + ggplot2::annotate(
            geom = "segment",
            x = rep(0,nrow(loading.supp)),
            xend = loading.supp[,comp[1]],
            y = rep(0,nrow(loading.supp)),
            yend = loading.supp[,comp[2]],
            col = "black",
            arrow = ggplot2::arrow(length = grid::unit(0.2, "cm")),size = 0.70
          ) + ggplot2::annotate(
            geom = "label",
            x = loading.supp[,comp[1]],
            y = loading.supp[,comp[2]],
            label = rownames(loading.supp),
            col = "blue",size = size.label
          )

      }

      return(graph.var)
    }

    #INDIVIDUAL
    if (choice == "ind"){
      limit.x = c(min(components[,comp[1]]), max(components[,comp[1]]))
      limit.y = c(min(components[,comp[2]]), max(components[,comp[2]]))

      Graph.observations <-
        ggplot2::ggplot(
          data.frame(components),
          ggplot2::aes(
            x = components[,comp[1]],
            y = components[,comp[2]],
            label = rownames(components)
          )
        ) +
        ggplot2::ggtitle("Factorial representation of individuals") +
        ggplot2::geom_point() +
        ggplot2::geom_text(ggplot2::aes(label=rownames(components), color = coloring.indiv),hjust=1, vjust=1,size = size.label)+
        ggplot2::xlab(paste(paste(nom.comp[1], inertie[comp[1],1]," %"))) +
        ggplot2::ylab(paste(paste(nom.comp[2], inertie[comp[2],1]," %"))) +
        ggplot2::scale_x_continuous(limits = as.numeric(limit.x)) +
        ggplot2::scale_y_continuous(limits = as.numeric(limit.y)) +
        ggplot2::geom_hline(
          yintercept = 0,
          linetype = "dotted",
          color = "black",
          size = 1
        ) +
        ggplot2::geom_vline(
          xintercept = 0,
          linetype = "dotted",
          color = "black",
          size = 1
        ) +
        ggplot2::theme_classic(base_size = size.legend)

      # if(choice == "qualitative" & supp.var == TRUE) {message("Modalities are represented as barycenter")}

      #Si il y a une variable qualitative supplementaire
      if (supp.var == TRUE){
        barycentre <- do.call(rbind.data.frame, res.PCAOS$Supp.var$barycenters)

        Graph.observations <-
          Graph.observations  + ggplot2::annotate(geom = "label",x = barycentre[,comp[1]], y = barycentre[,comp[2]],label =  rownames(barycentre),size = size.label,col = "blue")

        if(ellipse == TRUE){
          mat = list(NULL)

          for (var in which(level.scale.supp == 'nom' | level.scale.supp == 'ord')){

            var.quali = var.supp[,var]
            coord.ellipse <- cbind.data.frame(ellipses.coord(var.quali = var.quali,components = components[,comp],level.conf = level.conf))
            modalites <- levels(as.factor(var.quali))
            nb.modal <- length(modalites)

            x <- matrix(NA,nrow(coord.ellipse),nb.modal)
            y <- matrix(NA,nrow(coord.ellipse),nb.modal)

            compteur <- 1
            for (modal in seq(from = 1,
                              to = (nb.modal * 2),
                              by = 2)) {
              x[, compteur] <- coord.ellipse[, modal]
              y[, compteur] <- coord.ellipse[, modal + 1]
              compteur <- compteur + 1
            }

            mat[[var]] <- cbind(c(x),c(y))
            mat[[var]] <- cbind.data.frame(mat[[var]], as.vector(unlist(sapply(1:nb.modal, function(j) {rep(modalites[j],nrow(coord.ellipse))}))))


          }
          mat.d <- do.call(rbind.data.frame,mat)
          Graph.observations <- Graph.observations +  ggplot2::annotate(
            geom = "path",
            x =  mat.d[,1],
            y = mat.d[,2],
            size = 0.5,
            color = "darkblue",
            group = as.factor(mat.d[,3])
          )

        }

      }

      return(Graph.observations)
    }

    #QUANTIFICATION
    if (choice == "quantif"){
      graphs.list <- list(NA)
      modalities <- NULL
      for (var in 1:nb.var){
        if(level.scale[var] == "nom"){
          data.quantif <- data.frame(cbind(as.numeric(quantification[[var]]),rownames(quantification[[var]])))
          colnames(data.quantif) <- c("quantification","modalities")
          graphs.list[[var]] <- ggplot2::ggplot(data = data.quantif) +
            ggplot2::geom_point(ggplot2::aes(x=modalities, y=as.numeric(quantification)),size = 3) +
            #ggplot2::scale_y_continuous(limits=limit) +
            ggplot2::scale_x_discrete(limits = levels(data[,var])) +
            ggplot2::xlab("Categories") +
            ggplot2::ylab("Quantifications") +
            ggplot2::ggtitle(variables[var])+
            ggplot2::theme_classic(base_size = size.legend)
        }
        if (level.scale[var] == "ord"){
          data.quantif <- cbind(as.numeric(quantification[[var]]),rownames(quantification[[var]]))
          colnames(data.quantif) <- c("quantification","modalities")
          graphs.list[[var]] <- ggplot2::ggplot(data = data.frame(data.quantif), ggplot2::aes(x=modalities, y=as.numeric(quantification))) +
            ggplot2::geom_point(size = 3) +
            #ggplot2::scale_y_continuous(limits=limit) +
            ggplot2::scale_x_discrete(limits= levels(data[,var])) +
            ggplot2::xlab("Categories") +
            ggplot2::ylab("Quantifications") +
            ggplot2::ggtitle(variables[var]) +
            ggplot2::theme_classic(base_size = size.legend)
        }
        if (level.scale[var] == "num"){
          colnames(quantification[[var]]) <- c("values","quantification")
          quantification.num <- data.frame(quantification[[var]])

          graphs.list[[var]] <-
            ggplot2::ggplot(data = data.frame(quantification.num), ggplot2::aes(x = as.numeric(quantification.num[,1]),
                                                                                   y = as.numeric(quantification.num[,2]))) +
            ggplot2::geom_point(size = 3) +
            #ggplot2::scale_y_continuous(limits=limit) +
            ggplot2::xlab("Values") +
            ggplot2::ylab("Quantifications") +
            ggplot2::ggtitle(paste(variables[var]))+
            ggplot2::geom_line()+
            ggplot2::theme_classic(base_size = size.legend)
        }
      }

      if(is.null(sub.var.quantif)){
        plot.tot <- ggpubr::ggarrange(plotlist=graphs.list)
        return(ggpubr::annotate_figure(plot.tot,top = ggpubr::text_grob("Quantification of ..", color = "Black", face = "bold", size = 10)))

      }

      if(!is.null(sub.var.quantif)){
        if (length(sub.var.quantif) == 1){
          return(graphs.list[[sub.var.quantif]])
        }
        else{
          sub.list.graph <- list(NULL)
          c <- 1
          for (i in sub.var.quantif){
            sub.list.graph[[c]] <- graphs.list[[i]]
            c <- c + 1
          }
          return(ggpubr::ggarrange(plotlist=sub.list.graph))
        }
      }

    }

    #COMMON ELEMENTS
    if (choice == "qualitative" | choice == "all.var"){
      category.coord <- list(NULL)
      compteur <- 1
      #Calcul coordonnees des modalites (avec rk1)
      for (var in var.quali){
        category.coord[[compteur]] <- rep(NA,3)
        if(label.cat == 'var+cat'){
          labels.of.cat <- paste(colnames(data[,var,drop = FALSE]),rownames(quantification[[var]]),sep = "_")
        }
        if(label.cat == 'cat'){
          labels.of.cat <- paste(rownames(quantification[[var]]),sep = "_")
        }
        category.coord[[compteur]] <-
          cbind(
            labels.of.cat,
            as.numeric(quantification[[var]]) * weights[[var]][comp[1]],
            as.numeric(quantification[[var]]) * weights[[var]][comp[2]]
          )
        colnames(category.coord[[compteur]]) <- c("Modalites",nom.comp)
        category.coord[[compteur]] <- data.frame(category.coord[[compteur]])
        compteur <- compteur + 1
      }
      names(category.coord) <- colnames(data[,var.quali,drop = F])
      data.modal <- do.call(rbind.data.frame,category.coord)
      colnames(data.modal) <- colnames(category.coord[[1]])

      # #Calcul coordonnees des modalites (sans rk1)
      # if (res.PCAOS$Quantification$summary$rank == "no.restriction"){
      #   data.nom <- data.frame(data[,var.nom,drop = F])
      #   nom.multiple <- which(sapply(1:ncol(data.nom),function(x){length(table(data.nom[,x])) > 2}))
      #   data.nom.multiple <- data.nom[,nom.multiple,drop = F]
      #   modalite.multiple <- sapply(1:ncol(data.nom.multiple),function(var){levels(data.nom.multiple[,var])},simplify = FALSE)
      #   nb.modal.mulitple <- unlist(lapply(modalite.multiple,length))
      #
      #   for (var in 1:length(data.nom.multiple)){
      #     for (modal in 1:nb.modal[var]){
      #       labels.of.cat <- paste(colnames(data.nom.multiple[,var,drop = FALSE]),levels(as.factor(data.nom.multiple[,var,drop = FALSE])),sep = "_")
      #       category.coord[[nom.multiple[var]]] <- data.frame()
      #       category.coord[[nom.multiple[var]]][modal,c(2,3)] <- colMeans(components[which(data.nom.multiple[,var] == modalite.multiple[[var]][modal]),])[comp]
      #     }
      #   }
      # }

      #Construction of matrix mixed, containing all coordinates of all variables and categories
      #Qualitative var
      nb.modal <- unlist(lapply(category.coord,nrow))
      nb.var <- length(nb.modal)
      i.variables <- as.vector(unlist(sapply(1:nb.var, function(j) {rep(colnames(data[,var.quali,drop = F])[j],nb.modal[j])})))
      nature.quali <- NULL
      level.scale.quali <- level.scale[which(level.scale == 'nom'|level.scale == 'ord')]
      for (i in 1:length(nb.modal)){
        nature.quali <- c(nature.quali,rep(level.scale.quali[i],nb.modal[i]))
      }
      p.nom <- which(nature.quali == 'nom')
      p.ord <- which(nature.quali == 'ord')
      p.quali <-  1:nrow(data.modal)
      data.modal <- cbind(data.modal,i.variables,nature.quali)
      row.names(data.modal) <- data.modal[,1]

      #Calcul frequence de citation
      data.quali <- data.frame(data[,var.quali,drop = F])
      freq <- unlist(sapply(1:ncol(data.quali),function(var){table(data.quali[,var])},simplify = F))
      freq <- freq/nrow(data) * 100

      #Numeric var
      if(any(level.scale == 'num')){
        p.num <- 1:length(which(level.scale == "num"))
        p.nom <- p.nom + length(p.num)
        p.quali <-  p.quali + length(which(level.scale == "num"))
        weight.num <- data.frame(t(data.frame(weights)))
        weight.num <- weight.num[which(level.scale == "num"),]
        nb.var.num <- nrow(weight.num)
        weight.num <- cbind(weight.num[,comp],'num')
        colnames(weight.num) <- c(nom.comp[1],nom.comp[2],'nature')

        #Mixed var
        nature <- c(weight.num[,3],data.modal[,5])
        #na <- rep(NA,(sum(nb.modal) + nb.var))
        mixed <- data.frame()
        mixed <- rbind(weight.num[,c(1,2)],data.modal[,c(2,3)])
        mixed <- cbind(mixed,c(row.names(weight.num),i.variables),nature)
        mixed <- cbind(mixed,c(rep(NA,length(which(mixed[,4] == "num"))),freq))
      }else{
        nature <- c(data.modal[,5])
        mixed <- data.modal[,c(2,3)]
        mixed <- cbind(mixed,i,variables,nature)
        mixed <- cbind(mixed,freq)
        p.num <- 0
        nature <- c(data.modal[,5])
        #na <- rep(NA,(sum(nb.modal) + nb.var))
        mixed <- rbind(data.modal[,c(2,3)])
        mixed <- cbind(mixed,i.variables,nature)
      }

      legend <- 'right'
      if(label.size.freq == FALSE){
        mixed[,5] <- 12
        legend <- 'none'
      }

    }

    #MODALITIES REPRESENTATION
    if (choice == "qualitative"){

      nb.modal <- unlist(lapply(category.coord,nrow))
      nb.var <- length(nb.modal)
      identification.variable <- as.vector(unlist(sapply(1:nb.var, function(j) {rep(colnames(data[,var.quali])[j],nb.modal[j])})))
      data.modal <- data.frame(data.modal,identification.variable)
      data.modal <- data.frame(data.modal,freq)

      graph.modalites <- ggplot2::ggplot(data = data.modal)+
        ggplot2::geom_hline(
          yintercept = 0,
          linetype = "dotted",
          color = "black",
          size = 1
        ) +
        ggplot2::geom_vline(
          xintercept = 0,
          linetype = "dotted",
          color = "black",
          size = 1
        ) +
        ggplot2::ggtitle("Factorial representation of qualitative variables") +
        ggplot2::xlab(paste(nom.comp[1], inertie[comp[1], 1], " %")) +
        ggplot2::ylab(paste(nom.comp[2], inertie[comp[2], 1], " %")) +
        ggplot2::theme_classic(base_size = size.legend) +
        ggplot2::guides(size = ggplot2::guide_legend(title = "Citation frequency (%)"))


      if(ordinal.as.direction == TRUE & any(level.scale == 'ord')){
        data.ord <- mixed[which(mixed[,4] == 'ord'),]
        ord.var <- unique(data.ord[,3])
        x.max <- x.min <- rep(NA,length(ord.var))
        y.max <- y.min <- rep(NA,length(ord.var))
        for (i in 1:length(ord.var)){
          data.var <- data.ord[which(data.ord[,3] == ord.var[i]),]

          x.max[i] <- as.numeric(data.var[nrow(data.var),1])
          x.min[i] <- as.numeric(data.var[1,1])
          y.max[i] <- as.numeric(data.var[nrow(data.var),2])
          y.min[i] <- as.numeric(data.var[1,2])
        }

        graph.modalites <- graph.modalites +
          ggplot2::annotate(
            geom = "segment",
            x = x.min,
            xend = x.max,
            y = y.min,
            yend = y.max,
            col = "black",
            arrow = ggplot2::arrow(length = grid::unit(0.45, "cm")),size = 0.70
          )
      }

      if (label.size.freq == TRUE){
        graph.modalites <- graph.modalites +
          ggplot2::geom_label(
            ggplot2::aes(
              x = as.numeric(data.modal[, 2]),
              y =  as.numeric(data.modal[, 3]),
              label = data.modal[, 1],
              size = as.numeric(data.modal[,7])
            )
            #size = size.label
          )
      }else{
        graph.modalites <- graph.modalites +
          ggplot2::geom_label(
            ggplot2::aes(
              x = as.numeric(data.modal[, 2]),
              y =  as.numeric(data.modal[, 3]),
              label = data.modal[, 1]
            ),
            size = size.label
          )
      }
      if(supp.var == TRUE){
        if(any(level.scale.supp == "nom" ) | any(level.scale.supp == "ord" )){
          coord.supp.quali <- do.call(rbind.data.frame,res.PCAOS$Supp.var$coord.supp.quali)
          graph.modalites <-
            graph.modalites  + ggplot2::annotate(geom = "label",x = coord.supp.quali[,comp[1]], y = coord.supp.quali[,comp[2]],label =  rownames(coord.supp.quali),size = size.label,col = "blue")
        }
      }


      return(graph.modalites)

    }

    #MIXED VARIABLES
    if (choice == "all.var"){

      #invisible.var
      if (min.contribution > 0){
        weights <- t(data.frame(weights))
        invisible.var <- intersect(which(abs(as.numeric(weights[,1])) < min.contribution), which(abs(as.numeric(weights[,2])) < min.contribution))
        invisible.var <- rownames(weights[invisible.var,])
        row.to.supp <- which(mixed[,3] %in% invisible.var)
        if(length(row.to.supp) > 0){
          mixed <- mixed[-row.to.supp,]
          p.num <- which(mixed$nature == 'num')
          p.quali <- which(mixed$nature == 'quali')
        }
        print(paste('Variable',invisible.var,'is  not ploted'))
      }

      mix.graph <- ggplot2::ggplot() +
        ggplot2::geom_label(ggplot2::aes(
          x = as.numeric(mixed[p.num,1]),
          y = as.numeric(mixed[p.num,2]),
          label = rownames(mixed[p.num,])
        ),color = "black",size = size.label) +
        ggplot2::annotate(
          geom = "segment",
          x = rep(0,nrow(mixed[p.num,])),
          xend = as.numeric(mixed[p.num,1]),
          y = rep(0,nrow(mixed[p.num,])),
          yend = as.numeric(mixed[p.num,2]),
          col = "black",
          arrow = ggplot2::arrow(length = grid::unit(0.2, "cm")),size = 0.70
        ) +
        ggplot2::geom_hline(
          yintercept = 0,
          linetype = "dotted",
          color = "black",
          size = 1
        ) +
        ggplot2::geom_vline(
          xintercept = 0,
          linetype = "dotted",
          color = "black",
          size = 1
        ) +
        ggplot2::ggtitle("Factorial representation of all variables") +
        ggplot2::xlab(paste(nom.comp[1], inertie[comp[1],1]," %")) +
        ggplot2::ylab(paste(nom.comp[2], inertie[comp[2],1]," %")) +
        ggplot2::theme_classic(base_size = size.legend)  + ggplot2::guides(size = ggplot2::guide_legend(title = "Citation frequency (%)")) +
        ggplot2::theme(legend.position = legend)
      #+ scale_y_continuous(limits = c(min(mixed[,c(1,2)]),max(mixed[,c(1,2)])))

      if(ordinal.as.direction == TRUE & any(level.scale == 'ord')){
        data.ord <- mixed[which(mixed[,4] == 'ord'),]
        ord.var <- unique(data.ord[,3])
        x.max <- x.min <- rep(NA,length(ord.var))
        y.max <- y.min <- rep(NA,length(ord.var))
        for (i in 1:length(ord.var)){
          data.var <- data.ord[which(data.ord[,3] == ord.var[i]),]
          x.max[i] <- as.numeric(data.var[nrow(data.var),1])
          x.min[i] <- as.numeric(data.var[1,1])
          y.max[i] <- as.numeric(data.var[nrow(data.var),2])
          y.min[i] <- as.numeric(data.var[1,2])
        }

        mix.graph <- mix.graph +
          ggplot2::annotate(
            geom = "segment",
            x = x.min,
            xend = x.max,
            y = y.min,
            yend = y.max,
            col = "black",
            arrow = ggplot2::arrow(length = grid::unit(0.45, "cm")),size = 0.70
          )
      }

      if (label.size.freq == TRUE){
        mix.graph <- mix.graph +
          ggplot2::geom_text(ggplot2::aes(
            x = as.numeric(mixed[p.quali,1]),
            y =  as.numeric(mixed[p.quali,2]),
            label = rownames(mixed[p.quali,]),
            size = mixed[p.quali, 5]
          ),color = "black")
      }else{
        mix.graph <- mix.graph +
          ggplot2::geom_text(ggplot2::aes(
            x = as.numeric(mixed[p.quali,1]),
            y =  as.numeric(mixed[p.quali,2]),
            label = rownames(mixed[p.quali,])),
            size = size.label,color = "black")
      }

      if(supp.var == TRUE){
        if(any(level.scale.supp == "num" )){
          #Numeric variables
          loading.supp <- do.call(rbind.data.frame,coord.supp.num)
          mix.graph <-
            mix.graph  + ggplot2::annotate(
              geom = "segment",
              x = rep(0,nrow(loading.supp)),
              xend = loading.supp[,comp[1]],
              y = rep(0,nrow(loading.supp)),
              yend = loading.supp[,comp[2]],
              col = "black",
              arrow = ggplot2::arrow(length = grid::unit(0.2, "cm")),size = 0.70
            ) + ggplot2::annotate(
              geom = "label",
              x = loading.supp[,comp[1]],
              y = loading.supp[,comp[2]],
              label = rownames(loading.supp),
              col = "blue",size = size.label
            )
        }

        if(any(level.scale.supp == "nom" ) | any(level.scale.supp == "ord" )){
          coord.supp.quali <- do.call(rbind.data.frame,res.PCAOS$Supp.var$coord.supp.quali)
          mix.graph <-
            mix.graph  + ggplot2::annotate(geom = "label",x = coord.supp.quali[,comp[1]], y = coord.supp.quali[,comp[2]],label =  rownames(coord.supp.quali),size = size.label,col = "blue")
        }

      }

      return(mix.graph)
    }

  }
