#' PCAOS.plot
#'
#' Visualisation of results from PCAOS method. See details for available plot.
#'
#' @param res.PCAOS an object of class PCAOS
#' @param choice the graph to plot possible values are "screeplot","quantif","indiv","cor","modalities","mixed","squared loadings". See Details.
#' @param comp a length 2 vector with the components to plot
#' @param coloring.indiv A vector of length N to color individuals. If NULL, no coloring is applied.
#' @param sub.var.quantif a vector with variable of interest
#' @param supp.var TRUE or FALSE; if TRUE supplementary variables are added in factorial representation
#' @param conf.ellipsises boolean (FALSE by default), if TRUE, draw ellipses around categories of the qualitative variable as supplementary.
#' @param level.conf level of confidence ellipses
#' @param size.label size of label in graphs
#' @param size.legend size of label in graphs

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
#'   \item  qualitative: factorial representation of qualitatives variables trough the representation of it's categories. Coordinates of each category is calculted such as the single quantification of the category multiplied by the loading of the associated variable.
#' }
#' For mixed  variables
#'  \itemize{
#'   \item  mixed: factorial representation of mixed variables (weight for numeric variables, and categories for qualitative variables)
#'   \item  squared.loadings: plot of the squared loadings of all the variables.
#' }
#'
#' All graph are ggplot object
#'
#' @examples
#' data (antibiotic)
#' res.mix <- PCAOS::mix.data(antibiotic)
#' nature <- rep(NA,ncol(antibiotic)) #Setting nature argument
#' nature[res.mix$p.numeric] <- "num"
#' nature[res.mix$p.quali] <- "nom"
#' nature[c(1,15)] <- "ord"
#' nature
#'
#'res.PCAOS <- PCAOS::PCAOS(
#'  data = antibiotic,
#'  nature = nature,
#'  supp.var = c(1,2)
#')
#'
#'#Individuals graph
#'PCAOS::plot.PCAOS(res.PCAOS = res.PCAOS,choice = "ind",coloring = antibiotic$Atb.conso)
#'
#'#Quantifications
#'PCAOS::plot.PCAOS(res.PCAOS = res.PCAOS,choice = "quantif",sub.var.quantif =  c(4,8))
#'
#'#Qualitative variables
#'PCAOS::plot.PCAOS(res.PCAOS = res.PCAOS,choice = "qualitative",supp.var = TRUE,conf.ellipsises = TRUE)
#'
#'#Numeric variables
#'PCAOS::plot.PCAOS(res.PCAOS = res.PCAOS,choice = "numeric",supp.var = TRUE)
#'
#'#Mixed variables
#'PCAOS::plot.PCAOS(res.PCAOS = res.PCAOS,choice = "mixed",supp.var = TRUE,conf.ellipsises = TRUE)
#'
#' @export plot.PCAOS
#' @export
plot.PCAOS <-
  function(res.PCAOS,
           choice = "screeplot",
           comp = c(1,2),
           coloring.indiv = NULL,
           supp.var = FALSE,
           sub.var.quantif = NULL,
           conf.ellipsises = FALSE,
           level.conf = 0.95,
           size.label = 3.5,
           size.legend = 10
           ) {
    nb.comp <- ncol(res.PCAOS$components)
    check.plot.arg(choice,res.PCAOS$nature,res.PCAOS$nature.supp,supp.var,comp,nb.comp)

    #INFORMATION ABOUT THE RESULTS
    data <- res.PCAOS$data
    nature <- res.PCAOS$nature
    nom.comp <- c(paste("PC",comp[1],sep=""),paste("PC",comp[2],sep=""))
    inertie <- res.PCAOS$inertia

    #QUANTIFICATIONS
    graphs.list <- list(NA)
    variables <- colnames(data)
    nb.var <- length(variables)
    quantification <- vector("list", length = nb.var)
    quantification.nom <- res.PCAOS$quantification.categories.nom
    quantification.ord <- res.PCAOS$quantification.categories.ord
    var.nom <- names(quantification.nom)
    var.ord <- names(quantification.ord)
    var.num <- names(data[which(nature == "num")])
    compteur <- 1
    for (i in which(nature == "nom")) {
      quantification[[i]] <- res.PCAOS$quantification.categories.nom[[compteur]]
      compteur <- compteur + 1
    }
    compteur <- 1
    for (i in which(nature == "ord")) {
      quantification[[i]] <- res.PCAOS$quantification.categories.ord[[compteur]]
      compteur <- compteur + 1
    }
    compteur <- 1
    for (i in which(nature == "num")){
      quantification[[i]] <- cbind(data[,i],res.PCAOS$quantified.data[,i])
      compteur <- compteur + 1
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
    if (choice == "numeric"){
      data.graph.var <- data.frame(t(data.frame(res.PCAOS$weights)))
      weight.num <- data.graph.var[which(nature == "num"),]

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
    if (choice == "ind"){
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

      #Si il y a une variable qualitative supplementaire
      if (supp.var == TRUE){
        barycentre <-
          do.call(rbind.data.frame, res.PCAOS$coord.supp.quali)

        Graph.observations <-
          Graph.observations  + ggplot2::annotate(geom = "label",x = barycentre[,comp[1]], y = barycentre[,comp[2]],label =  rownames(barycentre),size = size.label,col = "blue")

        if(conf.ellipsises == TRUE){
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
      #LIMITS OF GRAPHS
      limit = c(min(res.PCAOS$quantified.data)-0.1, max(res.PCAOS$quantified.data)+0.1)
      for (var in 1:nb.var){
        if(nature[var] == "nom"){
          data.pour.graph <- cbind(as.numeric(quantification[[var]]),rownames(quantification[[var]]))
          colnames(data.pour.graph) <- c("quantification","modalities")
          graphs.list[[var]] <- ggplot2::ggplot(data = data.frame(data.pour.graph),
                                                ggplot2::aes(x=modalities, y=as.numeric(quantification))) +
            ggplot2::geom_point(size = 3) +
            ggplot2::scale_y_continuous(limits=limit) +
            ggplot2::xlab("Categories") +
            ggplot2::ylab("Quantifications") +
            ggplot2::ggtitle(variables[var])+
            ggplot2::theme_classic(base_size = size.legend)
        }
        if (nature[var] == "ord"){
          data.pour.graph <- cbind(as.numeric(quantification[[var]]),rownames(quantification[[var]]))
          colnames(data.pour.graph) <- c("quantification","modalities")
          graphs.list[[var]] <- ggplot2::ggplot(data = data.frame(data.pour.graph), ggplot2::aes(x=modalities, y=as.numeric(quantification))) +
            ggplot2::geom_point(size = 3) +
            ggplot2::scale_y_continuous(limits=limit) +
            ggplot2::xlab("Categories") +
            ggplot2::ylab("Quantifications") +
            ggplot2::ggtitle(variables[var]) +
            ggplot2::theme_classic(base_size = size.legend)
        }
        if (nature[var] == "num"){
          colnames(quantification[[var]]) <- c("quantification","values")
          graphs.list[[var]] <-
            ggplot2::ggplot(data = data.frame(quantification[[var]]), ggplot2::aes(x = as.numeric(quantification),
                                                                                   y = as.numeric(values))) +
            ggplot2::geom_point(size = 3) +
            ggplot2::scale_y_continuous(limits=limit) +
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
    if (choice =="qualitative" | choice == "mixed" ){
      category.coord <- list(NULL)
      compteur <- 1
      var.quali <- which(nature == "nom" | nature =="ord")
      nb.var.quali <- length(var.quali)

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

      max.x <- max(c(as.numeric(unlist(lapply(1:nb.var.quali, function(v)category.coord[[v]][,2]))),res.PCAOS$weights[[var]][1]))
      min.x <- min(c(as.numeric(unlist(lapply(1:nb.var.quali, function(v)category.coord[[v]][,2]))),res.PCAOS$weights[[var]][1]))

      max.y <- max(c(as.numeric(unlist(lapply(1:nb.var.quali, function(v)category.coord[[v]][,3]))),res.PCAOS$weights[[var]][2]))
      min.y <- min(c(as.numeric(unlist(lapply(1:nb.var.quali, function(v)category.coord[[v]][,3]))),res.PCAOS$weights[[var]][2]))

      limit.x <- c(min.x-0.2,max.x+0.2)
      limit.y <- c(min.y-0.2,max.y+0.2)

      data.modal <- category.coord[[1]]
      for (i in 2:length(category.coord)){
        data.modal <- rbind(data.modal,category.coord[[i]])
      }
      colnames(data.modal) <- colnames(category.coord[[1]])

    }

    #MODALITIES REPRESENTATION
    if (choice == "qualitative" ){
      nb.modal <- unlist(lapply(category.coord,nrow))
      nb.var <- length(nb.modal)
      identification.variable <- as.vector(unlist(sapply(1:nb.var, function(j) {rep(colnames(data[,var.quali])[j],nb.modal[j])})))
      data.modal <- data.frame(data.modal,identification.variable)

      graph.modalites <- ggplot2::ggplot(data = data.modal) +
        ggplot2::geom_label(ggplot2::aes(
          x = as.numeric(data.modal[,2]),
          y =  as.numeric(data.modal[,3]),
          label = data.modal[, 1],
          fill = as.character(identification.variable)
        ),
        size = size.label) +
        ggplot2::geom_line(ggplot2::aes(
          x = as.numeric(data.modal[,2]),
          y = as.numeric(data.modal[,3]),
          group = identification.variable,
          col = as.character(identification.variable)
        ),
        size = 1,show.legend = F) +
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
        ggplot2::theme_classic(base_size = size.legend) + ggplot2::guides(fill = ggplot2::guide_legend(title = "Variables",
                                                                                              override.aes = ggplot2::aes(label = "")))

      return(graph.modalites)

    }

    #MIXED VARIABLES
    if (choice == "mixed"){
      #NUM DATA
      data.graph.var <- data.frame(t(data.frame(res.PCAOS$weights)))
      weight.num <- data.graph.var[which(nature == "num"),]

      #NOM ET ORD DATA
      #Identification of the number of modalities per variable (for fill argument)
      nb.modal <- unlist(lapply(category.coord,nrow))
      nb.var <- length(nb.modal)
      identification.variable <- as.vector(unlist(sapply(1:nb.var, function(j) {rep(colnames(data[,var.quali])[j],nb.modal[j])})))
      data.modal <- data.frame(data.modal,identification.variable)


      mix.graph <- ggplot2::ggplot() +
        ggplot2::geom_label(ggplot2::aes(
          x = as.numeric(weight.num[,comp[1]]),
          y = as.numeric(weight.num[,comp[2]]),
          label = rownames(weight.num)
        ),color = "black",size = size.label) +
        ggplot2::annotate(
          geom = "segment",
          x = rep(0,nrow(weight.num)),
          xend = weight.num[,comp[1]],
          y = rep(0,nrow(weight.num)),
          yend = weight.num[,comp[2]],
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
        ggplot2::geom_label(ggplot2::aes(
          x = as.numeric(data.modal[,2]),
          y =  as.numeric(data.modal[,3]),
          label = data.modal[, 1],
          fill = as.character(identification.variable)
        ),
        size = size.label) +
        ggplot2::geom_line(ggplot2::aes(
          x = as.numeric(data.modal[,2]),
          y = as.numeric(data.modal[,3]),
          group = as.character(identification.variable),
          col = as.character(identification.variable)
        ),
        size = 1,show.legend = F) +
        ggplot2::ggtitle("Factorial representation of mixed variables") +
        ggplot2::xlab(paste(nom.comp[1], inertie[comp[1],1]," %")) +
        ggplot2::ylab(paste(nom.comp[2], inertie[comp[2],1]," %")) +
        ggplot2::theme_classic(base_size = size.legend)  + ggplot2::guides(fill = ggplot2::guide_legend(title = "Variables",
                                                                    override.aes = ggplot2::aes(label = "")))
      if(supp.var == TRUE){
        if(any(res.PCAOS$nature.supp == "num" )){
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
      }

      return(mix.graph)
    }


  }
