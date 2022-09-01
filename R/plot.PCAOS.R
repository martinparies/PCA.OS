#' plot.PCAOS
#'
#' Visualisation of results from PCAOS method. See details for available plot.
#'
#' @param res.PCAOS an object of class PCAOS
#' @param choice the available graphs are "screeplot","quantif","ind","numeric","qualitative","all.var". See Details.
#' @param comp a length 2 vector with the components to plot
#' @param coloring.indiv A vector of length N to color individuals. If NULL, no coloring is applied.
#' @param sub.var.quantif a vector with variable of interest
#' @param supp.var TRUE or FALSE; if TRUE supplementary variables are added in factorial representation
#' @param ellipse boolean (FALSE by default), if TRUE, draw ellipses around categories of the qualitative variable considered as supplementary.
#' @param level.conf level of confidence ellipses
#' @param size.label size of label in graphs
#' @param size.legend size of label in graphs
#' @param min.contribution variables with a contribution (i.e loading) lower than this value will not be plotted in the 'all.var' graph (useful for dataset with a lot of variables)

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
#'PCA.OS::plot.PCAOS(res.PCAOS = res.PCAOS,choice = "ind",coloring = antibiotic$Atb.conso)
#'PCA.OS::plot.PCAOS(res.PCAOS = res.PCAOS,choice = "ind",supp.var = TRUE,ellipse = TRUE)
#'
#'#Quantifications
#'PCA.OS::plot.PCAOS(res.PCAOS = res.PCAOS,choice = "quantif",sub.var.quantif =  c(4,8))
#'
#'#'#Numeric variables
#'PCA.OS::plot.PCAOS(res.PCAOS = res.PCAOS,choice = "qualitative",supp.var = TRUE)
#'
#'#Numeric variables
#'PCA.OS::plot.PCAOS(res.PCAOS = res.PCAOS,choice = "numeric",supp.var = TRUE)
#'
#'#All variables
#'PCA.OS::plot.PCAOS(res.PCAOS = res.PCAOS,choice = "all.var",supp.var = TRUE)
#'
#' @author
#' \itemize{
#'   \item Martin PARIES (Maintainer: \email{martin.paries@oniris-nantes.fr})
#'   \item Evelyne Vigneau
#'   \item Stephanie Bougeard
#' }
#'
#' @export plot.PCAOS
#' @export
#'
plot.PCAOS <-
  function(res.PCAOS,
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
           label.size.freq = TRUE
           ) {

    if(!inherits(res.PCAOS,"PCAOS")) stop("Non convenient object")

    nb.comp <- ncol(res.PCAOS$components)
    check.plot.arg(choice,res.PCAOS$level.scale,res.PCAOS$level.scale.supp,supp.var,comp,nb.comp,rank = res.PCAOS$summary$rank)

    #INFORMATION ABOUT THE RESULTS
    data <- res.PCAOS$data
    level.scale <- res.PCAOS$level.scale
    nom.comp <- c(paste("PC",comp[1],sep=""),paste("PC",comp[2],sep=""))
    inertie <- res.PCAOS$inertia

    #QUANTIFICATIONS
    variables <- colnames(data)
    nb.var <- length(variables)
    quantification <- vector("list", length = nb.var)
    quantification.nom <- res.PCAOS$quantification.categories.nom
    quantification.ord <- res.PCAOS$quantification.categories.ord
    var.nom <- names(quantification.nom)
    var.ord <- names(quantification.ord)
    var.num <- names(data[which(level.scale == "num")])
    compteur <- 1
    for (i in which(level.scale == "nom")) {
      quantification[[i]] <- res.PCAOS$quantification.categories.nom[[compteur]]
      compteur <- compteur + 1
    }
    compteur <- 1
    for (i in which(level.scale == "ord")) {
      quantification[[i]] <- res.PCAOS$quantification.categories.ord[[compteur]]
      compteur <- compteur + 1
    }
    compteur <- 1
    if(res.PCAOS$summary$rank == "one"){
      for (i in which(level.scale == "num")){
        quantification[[i]] <- cbind(data[,i],res.PCAOS$quantified.data[,i])
        compteur <- compteur + 1
      }
    }
    if(res.PCAOS$summary$rank == "no.restriction"){
      for (i in which(level.scale == "num")){
        quantification[[i]] <- cbind(data[,i],res.PCAOS$quantified.data[[i]])
        compteur <- compteur + 1
      }
    }
    names (quantification) <- variables

    #SCREEPLOT
    if(choice == "screeplot"){
      screeplot <-
        ggplot2::ggplot(inertie, ggplot2::aes(y = inertia, x = (1:nrow(inertie)))) +
        ggplot2::scale_x_discrete(labels = "") +
        ggplot2::geom_bar(stat = "identity", width = 0.5, fill = "steelblue") +
        ggplot2::geom_text(
          ggplot2::aes(label = paste(inertia, "%")),
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
    if (choice == "numeric" | (choice == "all.var" & (identical (level.scale,rep("num",length(level.scale)))))){
      data.graph.var <- data.frame(t(data.frame(res.PCAOS$weights)))
      weight.num <- data.graph.var[which(level.scale == "num"),]

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


      #Si il y a une variable qualitative supplementaire
      if (supp.var == TRUE){
        loading.supp <- do.call(rbind.data.frame,res.PCAOS$coord.supp.num)

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
    if (choice == "ind" | (choice == "qualitative" & (supp.var == TRUE & all(res.PCAOS$level.scale.supp == "nom")))){
      nom.indiv <- rownames(res.PCAOS$data)
      limit.x = c(min(res.PCAOS$components[,comp[1]]), max(res.PCAOS$components[,comp[1]]))
      limit.y = c(min(res.PCAOS$components[,comp[2]]), max(res.PCAOS$components[,comp[2]]))

      Graph.observations <-
        ggplot2::ggplot(
          data.frame(res.PCAOS$components),
          ggplot2::aes(
            x = res.PCAOS$components[,comp[1]],
            y = res.PCAOS$components[,comp[2]],
            label = rownames(res.PCAOS$components)
          )
        ) +
        ggplot2::ggtitle("Factorial representation of individuals") +
        ggplot2::geom_point() +
        ggplot2::geom_text(ggplot2::aes(label=rownames(res.PCAOS$components), color = coloring.indiv),hjust=1, vjust=1,size = size.label)+
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
        barycentre <-
          do.call(rbind.data.frame, res.PCAOS$coord.supp.quali)

        Graph.observations <-
          Graph.observations  + ggplot2::annotate(geom = "label",x = barycentre[,comp[1]], y = barycentre[,comp[2]],label =  rownames(barycentre),size = size.label,col = "blue")

        if(ellipse == TRUE){
          mat = list(NULL)

          for (var.supp in 1:ncol(res.PCAOS$quali.var.supp)){

            var.quali = res.PCAOS$quali.var.supp[,var.supp]
            coord.ellipse <- cbind.data.frame(ellipses.coord(var.quali = var.quali,components = res.PCAOS$components[,comp],level.conf = level.conf))
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

            mat[[var.supp]] <- cbind(c(x),c(y))
            mat[[var.supp]] <- cbind.data.frame(mat[[var.supp]], as.vector(unlist(sapply(1:nb.modal, function(j) {rep(modalites[j],nrow(coord.ellipse))}))))


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
      limit = c(min(res.PCAOS$quantified.data)-0.1, max(res.PCAOS$quantified.data)+0.1)
      for (var in 1:nb.var){
        if(level.scale[var] == "nom"){
          data.pour.graph <- cbind(as.numeric(quantification[[var]]),rownames(quantification[[var]]))
          colnames(data.pour.graph) <- c("quantification","modalities")
          graphs.list[[var]] <- ggplot2::ggplot(data = data.frame(data.pour.graph),
                                                ggplot2::aes(x=modalities, y=as.numeric(quantification))) +
            ggplot2::geom_point(size = 3) +
            #ggplot2::scale_y_continuous(limits=limit) +
            ggplot2::xlab("Categories") +
            ggplot2::ylab("Quantifications") +
            ggplot2::ggtitle(variables[var])+
            ggplot2::theme_classic(base_size = size.legend)
        }
        if (level.scale[var] == "ord"){
          data.pour.graph <- cbind(as.numeric(quantification[[var]]),rownames(quantification[[var]]))
          colnames(data.pour.graph) <- c("quantification","modalities")
          graphs.list[[var]] <- ggplot2::ggplot(data = data.frame(data.pour.graph), ggplot2::aes(x=modalities, y=as.numeric(quantification))) +
            ggplot2::geom_point(size = 3) +
            #ggplot2::scale_y_continuous(limits=limit) +
            ggplot2::xlab("Categories") +
            ggplot2::ylab("Quantifications") +
            ggplot2::ggtitle(variables[var]) +
            ggplot2::theme_classic(base_size = size.legend)
        }
        if (level.scale[var] == "num"){
          colnames(quantification[[var]]) <- c("quantification","values")
          graphs.list[[var]] <-
            ggplot2::ggplot(data = data.frame(quantification[[var]]), ggplot2::aes(x = as.numeric(quantification),
                                                                                   y = as.numeric(values))) +
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
    if (choice == "qualitative" | choice == "all.var" & any(level.scale == "nom" | level.scale == "ord")){
      category.coord <- list(NULL)
      compteur <- 1
      var.quali <- which(level.scale == "nom" | level.scale =="ord")
      nb.var.quali <- length(var.quali)

      if (res.PCAOS$summary$rank == "one"){
        for (var in var.quali){
          category.coord[[compteur]] <- rep(NA,3)
          category.coord[[compteur]] <-
            cbind(
              paste(colnames(data[,var,drop = FALSE]),rownames(quantification[[var]]),sep = "_"),
              as.numeric(quantification[[var]]) * res.PCAOS$weights[[var]][comp[1]],
              as.numeric(quantification[[var]]) * res.PCAOS$weights[[var]][comp[2]]
            )
          colnames(category.coord[[compteur]]) <- c("Modalites",nom.comp)
          category.coord[[compteur]] <- data.frame(category.coord[[compteur]])
          compteur <- compteur + 1
        }
        names(category.coord) <- colnames(data[,var.quali])

        data.modal <- category.coord[[1]]
        for (i in 2:length(category.coord)){
          data.modal <- rbind(data.modal,category.coord[[i]])
        }
        colnames(data.modal) <- colnames(category.coord[[1]])
      }

      if (res.PCAOS$summary$rank == "no.restriction"){
        var.nom <- which(level.scale == "nom" )
        data.quali <- data.frame(res.PCAOS$data[,var.nom,drop = F])
        for (i in 1:ncol(data.quali)){data.quali[,i] <- factor(data.quali[,i])}
        variables.quali <- colnames(data.quali)
        modalite <- sapply(1:ncol(data.quali),function(var){levels(data.quali[,var])},simplify = FALSE)
        nb.modal <- unlist(lapply(modalite,length))
        category.coord <- sapply(1:length(var.nom), function(var) matrix(NA,nb.modal[var],nb.comp),simplify = FALSE)

        for (var in 1:length(var.nom)){
          for (modal in 1:nb.modal[var]){
            category.coord[[var]][modal,] <- colMeans(res.PCAOS$components[which(data.quali[,var] == modalite[[var]][modal]),])
          }
          colnames(category.coord[[var]]) <- colnames(res.PCAOS$components)
          rownames(category.coord[[var]]) <- modalite[[var]]
          #category.coord[[var]] <- category.coord[[var]][order(rownames(category.coord[[var]])),]
        }
        names(category.coord) <- variables.quali
        # fill.arg <- NULL
        # for (var in 1:length(var.quali)){fill.arg <- c(fill.arg,rep (variables.quali[var],nb.modal[var]))}
        modalities <- unlist(sapply(1:length(var.nom),function(var){paste(variables.quali[var],modalite[[var]],sep = "_")},simplify = F))
        category.coord.tot <- do.call("rbind", category.coord)
        data.modal <- data.frame(modalities = modalities,category.coord.tot[,c(comp[1],comp[2])])

        var.ord <- which(level.scale == "ord" )
        data.quali <- data.frame(res.PCAOS$data[,var.ord,drop = F])

        category.coord.ord <- list(NULL)
        for (var in var.ord){
          category.coord.ord[[compteur]] <- rep(NA,3)
          category.coord.ord[[compteur]] <-
            cbind(
              paste(colnames(data[,var,drop = FALSE]),rownames(quantification[[var]]),sep = "_"),
              as.numeric(quantification[[var]]) * res.PCAOS$weights[[var]][comp[1]],
              as.numeric(quantification[[var]]) * res.PCAOS$weights[[var]][comp[2]]
            )
          colnames(category.coord.ord[[compteur]]) <- c("Modalites",nom.comp)
          category.coord.ord[[compteur]] <- data.frame(category.coord.ord[[compteur]])
          compteur <- compteur + 1
        }
        category.coord.ord <- do.call(rbind.data.frame,category.coord.ord)
        #rownames(category.coord.ord) <- category.coord.ord[,1]
        category.coord.ord <- category.coord.ord[,-1,drop = F]
        data.modal <- rbind(data.modal,category.coord.ord)
        colnames(data.modal) <- c("Modalites",nom.comp)
      }

      max.x <- as.numeric(max(data.modal[,2]))
      min.x <- as.numeric(min(data.modal[,2]))

      max.y <- as.numeric(max(data.modal[,3]))
      min.y <- as.numeric(min(data.modal[,3]))

      limit.x <- c(min.x-0.2,max.x+0.2)
      limit.y <- c(min.y-0.2,max.y+0.2)

      #Construction of matrix mixed, containing all coordinates of all variables and categories
      #quali var
      nb.modal <- unlist(lapply(category.coord,nrow))
      nb.var <- length(nb.modal)

      #Variables
      i.variables <- as.vector(unlist(sapply(1:nb.var, function(j) {rep(colnames(res.PCAOS$data[,var.quali,drop = F])[j],nb.modal[j])})))
      #data.modal <- cbind(data.modal,fill.arg)

      #Quali var
      p.nom <- 1:sum(nb.modal)
      data.modal <- cbind(data.modal,i.variables,'quali')
      row.names(data.modal) <- data.modal[,1]

      #Numeric var
      p.num <- 1:length(which(level.scale == "num"))
      p.nom <- p.nom + length(p.num)
      weight.num <- data.frame(t(data.frame(res.PCAOS$weights)))
      weight.num <- weight.num[which(level.scale == "num"),]
      weight.num <- cbind(weight.num[,comp],'num')
      nb.var.num <- nrow(weight.num)
      colnames(weight.num) <- c(nom.comp[1],nom.comp[2],'nature')

      #Mixed var
      nature <- c(weight.num[,3],data.modal[,5])
      na <- rep(NA,(sum(nb.modal) + nb.var))
      mixed <- data.frame()
      mixed <- rbind(weight.num[,c(1,2)],data.modal[,c(2,3)])
      mixed <- cbind(mixed,c(row.names(weight.num),i.variables),nature)

      #Calcul frequence de citation
      data.quali <- data.frame(res.PCAOS$data[,var.quali,drop = F])
      freq <- unlist(sapply(1:ncol(data.quali),function(var){table(data.quali[,var])},simplify = F))
      freq <- freq/nrow(res.PCAOS$data) * 100

      mixed <- cbind(mixed,c(rep(NA,length(which(mixed[,4] == "num"))),freq))

      legend <- 'right'
      if(label.size.freq == FALSE){
        mixed[,5] <- 1
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

      graph.modalites <- ggplot2::ggplot(data = data.modal) +
          ggplot2::geom_label(
          ggplot2::aes(
            x = as.numeric(data.modal[, 2]),
            y =  as.numeric(data.modal[, 3]),
            label = data.modal[, 1],
            size = as.numeric(data.modal[,5])
          ),
          #size = size.label
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
        ggplot2::ggtitle("Factorial representation of qualitative variables") +
        ggplot2::xlab(paste(nom.comp[1], inertie[comp[1], 1], " %")) +
        ggplot2::ylab(paste(nom.comp[2], inertie[comp[2], 1], " %")) +
        ggplot2::theme_classic(base_size = size.legend) +
        ggplot2::guides(size = ggplot2::guide_legend(title = "Citation frequency (%)"))

      # if(res.PCAOS$summary$rank == "one"){
      #   graph.modalites <- graph.modalites + ggplot2::geom_line(ggplot2::aes(
      #     x = as.numeric(data.modal[,2]),
      #     y = as.numeric(data.modal[,3]),
      #     group = identification.variable,
      #     col = as.character(identification.variable)
      #   ),
      #   size = 1,show.legend = F)
      # }

      if(supp.var == TRUE){
        if(any(res.PCAOS$level.scale.supp == "num" )){
          loading.supp <- do.call(rbind.data.frame,res.PCAOS$coord.supp.num)
          graph.modalites <-
            graph.modalites  + ggplot2::annotate(
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

        if(any(res.PCAOS$level.scale.supp == "nom" ) | any(res.PCAOS$level.scale.supp == "ord" )){
          barycentre <-
            do.call(rbind.data.frame, res.PCAOS$coord.supp.quali)

          graph.modalites <-
            graph.modalites  + ggplot2::annotate(geom = "label",x = barycentre[,comp[1]], y = barycentre[,comp[2]],label =  rownames(barycentre),size = size.label,col = "blue")

          if(ellipse == TRUE){
            mat = list(NULL)

            for (var.supp in 1:ncol(res.PCAOS$quali.var.supp)){

              var.quali = res.PCAOS$quali.var.supp[,var.supp]
              coord.ellipse <- cbind.data.frame(ellipses.coord(var.quali = var.quali,components = res.PCAOS$components,level.conf = level.conf))
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

              mat[[var.supp]] <- cbind(c(x),c(y))
              mat[[var.supp]] <- cbind.data.frame(mat[[var.supp]], as.vector(unlist(sapply(1:nb.modal, function(j) {rep(modalites[j],nrow(coord.ellipse))}))))


            }
            mat.d <- do.call(rbind.data.frame,mat)
            graph.modalites <- graph.modalites +  ggplot2::annotate(
              geom = "path",
              x =  mat.d[,1],
              y = mat.d[,2],
              size = 0.5,
              color = "darkblue",
              group = as.factor(mat.d[,3])
            )

          }
        }

      }



      return(graph.modalites)

    }

    #MIXED VARIABLES
    if (choice == "all.var"){

      # #NUM DATA
      # data.graph.var <- data.frame(t(data.frame(res.PCAOS$weights)))
      # weight.num <- data.graph.var[which(level.scale == "num"),]
      #
      # #NOM ET ORD DATA
      # #Identification of the number of modalities per variable (for fill argument)
      # nb.modal <- unlist(lapply(category.coord,nrow))
      # nb.var <- length(nb.modal)
      # identification.variable <- as.vector(unlist(sapply(1:nb.var, function(j) {rep(colnames(data[,var.quali])[j],nb.modal[j])})))
      # data.modal <- data.frame(data.modal,identification.variable)
      # data.modal <- data.frame(data.modal,freq)

      #invisible.var
      if (min.contribution > 0){
        weights <- t(data.frame(res.PCAOS$weights))
        invisible.var <- intersect(which(abs(as.numeric(weights[,1])) < min.contribution), which(abs(as.numeric(weights[,2])) < min.contribution))
        invisible.var <- rownames(weights[invisible.var,])
        row.to.supp <- which(mixed[,3] %in% invisible.var)
        if(length(row.to.supp) > 0){
          mixed <- mixed[-row.to.supp,]
          p.num <- which(mixed$nature == 'num')
          p.nom <- which(mixed$nature == 'quali')
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
        ggplot2::geom_text(ggplot2::aes(
          x = as.numeric(mixed[p.nom,1]),
          y =  as.numeric(mixed[p.nom,2]),
          label = rownames(mixed[p.nom,]),
          size = mixed[p.nom, 5]
        ),color = "black") +
        ggplot2::ggtitle("Factorial representation of all variables") +
        ggplot2::xlab(paste(nom.comp[1], inertie[comp[1],1]," %")) +
        ggplot2::ylab(paste(nom.comp[2], inertie[comp[2],1]," %")) +
        ggplot2::theme_classic(base_size = size.legend)  + ggplot2::guides(size = ggplot2::guide_legend(title = "Citation frequency (%)")) +
        ggplot2::theme(legend.position = legend)

      if(supp.var == TRUE){
        if(any(res.PCAOS$level.scale.supp == "num" )){
          #Numeric variables
          loading.supp <- do.call(rbind.data.frame,res.PCAOS$coord.supp.num)
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
        if(any(res.PCAOS$level.scale.supp == "nom" ) | any(res.PCAOS$level.scale.supp == "ord" )){
          barycentre <-
            do.call(rbind.data.frame, res.PCAOS$coord.supp.quali)

          mix.graph <-
            mix.graph  + ggplot2::annotate(geom = "label",x = barycentre[,comp[1]], y = barycentre[,comp[2]],label =  rownames(barycentre),size = size.label,col = "blue")

          if(ellipse == TRUE){
            mat = list(NULL)

            for (var.supp in 1:ncol(res.PCAOS$quali.var.supp)){

              var.quali = res.PCAOS$quali.var.supp[,var.supp]
              coord.ellipse <- cbind.data.frame(ellipses.coord(var.quali = var.quali,components = res.PCAOS$components,level.conf = level.conf))
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

              mat[[var.supp]] <- cbind(c(x),c(y))
              mat[[var.supp]] <- cbind.data.frame(mat[[var.supp]], as.vector(unlist(sapply(1:nb.modal, function(j) {rep(modalites[j],nrow(coord.ellipse))}))))


            }
            mat.d <- do.call(rbind.data.frame,mat)
            mix.graph <- mix.graph +  ggplot2::annotate(
              geom = "path",
              x =  mat.d[,1],
              y = mat.d[,2],
              size = 0.5,
              color = "darkblue",
              group = as.factor(mat.d[,3])
            )

          }
        }

      }

      return(mix.graph)
    }

  }
