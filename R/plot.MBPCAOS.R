#' plot.MBPCAOS
#'
#' Visualisation of results from MBPCAOS method. See details for available plot.
#'
#' @param res.MBPCAOS an object of class MBPCAOS
#' @param choice the graph to plot possible values are "screeplot","quantif","indiv","cor","modalities","mixed","squared loadings". See Details.
#' @param comp a length 2 vector with the components to plot
#' @param coloring.ind A vector of length N to color individuals. If NULL, no coloring is applied.
#' @param coloring.cat Categories are colored according to their rescective blocs (default = "blocs"). If "variables", categories are colored according to their variables.
#' @param sub.var.quantif a vector with variable of interest
#' @param supp.var boolean (FALSE by default), if TRUE supplementary variables are added in factorial representation

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
#' For blocs
#'  \itemize{
#'   \item  blocs: factorial representation of blocs, trough the correlation of block component with principal components
#' }
#'
#' All graph are ggplot object
#'
#' @export plot.MBPCAOS
#' @export
plot.MBPCAOS <-function(
  res.MBPCAOS,
  choice,
  comp = c(1,2),
  coloring.ind = NULL,
  coloring.cat = "blocs",
  supp.var = FALSE,
  sub.var.quantif = NULL) {

  nb.comp <- ncol(res.MBPCAOS$components)
  check.plot.MB.arg(choice,res.MBPCAOS$nature,res.MBPCAOS$nature.supp,supp.var,comp,nb.comp)


  #Information about the results
  nature <- res.MBPCAOS$nature
  variables <- colnames(res.MBPCAOS$data)
  nb.comp <- ncol(res.MBPCAOS$T.tot)
  nb.var <- length(variables)
  nb.bloc <- length(res.MBPCAOS$blocs.name)
  quantification <- vector("list", length = nb.var)
  quantification.nom <- res.MBPCAOS$quantification.categories.nom
  quantification.ord <- res.MBPCAOS$quantification.categories.ord
  var.nom <- names(quantification.nom)
  var.ord <- names(quantification.ord)
  var.num <- names(res.MBPCAOS$data[which(nature == "num")])
  inertie <- res.MBPCAOS$inertia
  limit = c(min(res.MBPCAOS$quantified.data), max(res.MBPCAOS$quantified.data))

  compteur <- 1
  for (i in which(nature == "nom")) {
    quantification[[i]] <- res.MBPCAOS$quantification.categories.nom[[compteur]]
    compteur <- compteur + 1
  }
  compteur <- 1
  for (i in which(nature == "ord")) {
    quantification[[i]] <- res.MBPCAOS$quantification.categories.ord[[compteur]]
    compteur <- compteur + 1
  }
  compteur <- 1
  for (i in which(nature == "num")){
    quantification[[i]] <- cbind(res.MBPCAOS$data[,i],res.MBPCAOS$quantified.data[,i])
    compteur <- compteur + 1
  }
  names (quantification) <- variables

  nom.comp <- c(paste("CP",comp[1],sep=""),paste("CP",comp[2],sep=""))

  #Screeplot
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
      ggplot2::theme_classic(base_size = 18)

    return (screeplot)
  }

  #Quantification
  if (choice == "quantif"){
    graphs.list <- list(NULL)
    for (var in 1:nb.var){
      if(nature[var] == "nom"){
        data.pour.graph <- cbind(as.numeric(quantification[[var]]),row.names(quantification[[var]]))
        colnames(data.pour.graph) <- c("quantification","modalities")
        graphs.list[[var]] <- ggplot2::ggplot(data = data.frame(data.pour.graph), ggplot2::aes(x=modalities, y=as.numeric(quantification))) +
          ggplot2::geom_point(size = 4) +
          ggplot2::scale_y_continuous(limits=limit) +
          ggplot2::xlab("Original data") +
          ggplot2::ylab("Quantifications") +
          ggplot2::ggtitle(variables[var])+
          ggplot2::theme_classic(base_size = 18)
      }
      if (nature[var] == "ord"){
        data.pour.graph <- cbind(as.numeric(quantification[[var]]),row.names(quantification[[var]]))
        colnames(data.pour.graph) <- c("quantification","modalities")
        graphs.list[[var]] <- ggplot2::ggplot(data = data.frame(data.pour.graph), ggplot2::aes(x=modalities, y=as.numeric(quantification))) +
          ggplot2::geom_point(size = 4) +
          #geom_line(aes(group = 1),linetype = "dashed", color = "black", size = 1) +
          ggplot2::scale_y_continuous(limits=limit) +
          ggplot2::xlab("Original data") +
          ggplot2::ylab("Quantifications") +
          ggplot2::ggtitle(variables[var]) +
          ggplot2::theme_classic(base_size = 18)
      }
      if (nature[var] == "num"){
        graphs.list[[var]] <- ggplot2::ggplot(data = data.frame(quantification[[var]]), ggplot2::aes(x=quantification[[var]][,1], y=quantification[[var]][,2])) +
          ggplot2::geom_point(size = 4) +
          ggplot2::scale_y_continuous(limits=limit) +
          ggplot2::xlab("Original data") +
          ggplot2::ylab("Quantifications") +
          ggplot2::ggtitle(paste(variables[var]))+
          ggplot2::geom_line()+
          ggplot2::theme_classic(base_size = 18)
      }
    }
    if(is.null(sub.var.quantif)){
      plot.tot <- ggpubr::ggarrange(plotlist=graphs.list)
      return(ggpubr::annotate_figure(plot.tot,top = ggpubr::text_grob("Quantification of ..", color = "Black", face = "bold", size = 15)))

    }

    if(!is.null(sub.var.quantif)){
      if (length(sub.var.quantif) == 1){
        print(graphs.list[[sub.var.quantif]])
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

  #Individuals
  if (choice == "ind"){
    nom.indiv <- rownames(res.MBPCAOS$JDD)
    limit.x = c(min(res.MBPCAOS$components[,comp[1]]), max(res.MBPCAOS$components[,comp[1]]))
    limit.y = c(min(res.MBPCAOS$components[,comp[2]]), max(res.MBPCAOS$components[,comp[2]]))
    Graph.observations <-
      ggplot2::ggplot(
        data.frame(res.MBPCAOS$components),
        ggplot2::aes(
          x = res.MBPCAOS$components[, comp[1]],
          y = res.MBPCAOS$components[, comp[2]],
          label = rownames(res.MBPCAOS$components)
        )
      ) +
      ggplot2::ggtitle("Factorial representation of individuals") +
      ggplot2::geom_point() +
      ggplot2::geom_text(ggplot2::aes(label = rownames(res.MBPCAOS$components), color = coloring.ind),
                hjust = 1,
                vjust = 1) +
      ggplot2::xlab(paste(nom.comp[1], res.MBPCAOS$inertie[comp[1], 2], " %")) +
      ggplot2::ylab(paste(nom.comp[2], res.MBPCAOS$inertie[comp[1], 2], " %")) +
      ggplot2::scale_x_continuous(limits = as.numeric(limit.x)) +
      ggplot2::scale_y_continuous(limits = as.numeric(limit.y)) +
      ggplot2::geom_hline(
        yintercept = 0,
        linetype = "dotted",
        color = "black",
        size = 1
      ) + ggplot2::geom_vline(
        xintercept = 0,
        linetype = "dotted",
        color = "black",
        size = 1
      ) + ggplot2::theme_classic(base_size = 18)
    return(Graph.observations)
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
          paste(colnames(res.MBPCAOS$data[,var,drop = FALSE]),rownames(quantification[[var]]),sep = "_"),
          as.numeric(quantification[[var]]) * res.MBPCAOS$weights[[var]][comp[1]],
          as.numeric(quantification[[var]]) * res.MBPCAOS$weights[[var]][comp[2]]
        )
      colnames(category.coord[[compteur]]) <- c("Modalites",nom.comp)
      category.coord[[compteur]] <- data.frame(category.coord[[compteur]])
      compteur <- compteur + 1
    }
    names(category.coord) <- colnames(res.MBPCAOS$data[,var.quali])

    max.x <- max(c(as.numeric(unlist(lapply(1:nb.var.quali, function(v)category.coord[[v]][,2]))),res.MBPCAOS$weights[[var]][1]))
    min.x <- min(c(as.numeric(unlist(lapply(1:nb.var.quali, function(v)category.coord[[v]][,2]))),res.MBPCAOS$weights[[var]][1]))

    max.y <- max(c(as.numeric(unlist(lapply(1:nb.var.quali, function(v)category.coord[[v]][,3]))),res.MBPCAOS$weights[[var]][2]))
    min.y <- min(c(as.numeric(unlist(lapply(1:nb.var.quali, function(v)category.coord[[v]][,3]))),res.MBPCAOS$weights[[var]][2]))

    limit.x <- c(min.x-0.2,max.x+0.2)
    limit.y <- c(min.y-0.2,max.y+0.2)

    data.modal <- category.coord[[1]]
    for (i in 2:length(category.coord)){
      data.modal <- rbind(data.modal,category.coord[[i]])
    }
    colnames(data.modal) <- colnames(category.coord[[1]])

  }

  #Qualitative variables
  if (choice == "qualitative" ){
    nb.modal <- unlist(lapply(category.coord,nrow))
    nb.var <- length(nb.modal)
    identification.variable <- as.vector(unlist(sapply(1:nb.var, function(j) {rep(colnames(res.MBPCAOS$data[,var.quali])[j],nb.modal[j])})))
    data.modal <- data.frame(data.modal,identification.variable)

    #Identification of the number of modalities per variable (for fill argument)
    if (coloring.cat == "variables"){
      fill.arg <- as.vector(unlist(sapply(1:nb.var, function(j) {rep(colnames(res.MBPCAOS$data[,var.quali])[j],nb.modal[j])})))
      data.modal <- cbind(data.modal,fill.arg)
      legend <- "Variables"
    }

    if (coloring.cat == "blocs"){
      identification.blocs <- unlist(sapply(1:nb.bloc,function(x)rep(res.MBPCAOS$blocs.name[x],res.MBPCAOS$blocs[x],)))
      identification.blocs <- identification.blocs[which(nature == "nom" | nature =="ord")]
      fill.arg <- unlist(sapply(1:nb.var,function(x)rep(identification.blocs[x],nb.modal[x],)))
      #variables.quali <- variables[which(nature == "nom" | nature =="ord")]
      #variables.modal <- unlist(sapply(1:nb.var,function(x)rep(variables.quali[x],nb.modal[x],)))
      #data.modal[,1] <- paste(variables.modal,data.modal[,1],sep="_")
      data.modal <- cbind(data.modal,fill.arg)
      legend <- "Blocs"
    }

    graph.modalites <- ggplot2::ggplot() +
      ggplot2::geom_label(ggplot2::aes(
        x = as.numeric(data.modal[,2]),
        y =  as.numeric(data.modal[,3]),
        label = data.modal[, 1],
        fill = as.character(fill.arg)
      ),
      size = 5) +
      ggplot2::geom_line(ggplot2::aes(
        x = as.numeric(data.modal[,2]),
        y = as.numeric(data.modal[,3]),
        group = identification.variable
      ),
      size = 0.5,show.legend = F) +
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
      ggplot2::xlab(paste(nom.comp[1], inertie[comp[1], 2], " %")) +
      ggplot2::ylab(paste(nom.comp[2], inertie[comp[2], 2], " %")) +
      ggplot2::theme_classic(base_size = 18) + ggplot2::guides(fill = ggplot2::guide_legend(title = "Variables",
                                                                                            override.aes = ggplot2::aes(label = "")))

    #Si il y a une variable qualitative supplementaire
    if (supp.var == TRUE){
      barycentre <-
        do.call(rbind.data.frame, res.MBPCAOS$coord.supp.quali)

      graph.modalites <-
        graph.modalites  + ggplot2::annotate(geom = "label",x = barycentre[,comp[1]], y = barycentre[,comp[2]],label =  rownames(barycentre),size = 5,col = "blue")

    }
    return(graph.modalites)
  }

  #Numeric variables
  if (choice == "numeric"){
    weight <- data.frame(t(data.frame(res.MBPCAOS$weights)))
    weight.num <- weight[which(nature == "num"),]
    identification.blocs <- unlist(sapply(1:nb.bloc,function(x)rep(res.MBPCAOS$blocs.name[x],res.MBPCAOS$blocs[x],)))
    identification.blocs <- identification.blocs[which(nature == "num")]

    facto.num <-
      ggplot2::ggplot(data = data.frame(weight.num),
                      ggplot2::aes (
                        x = as.numeric(weight.num[,comp[1]]),
                        y = as.numeric(weight.num[,comp[2]])
                      )) +
      ggplot2::geom_label(ggplot2::aes(label = row.names(weight.num),fill = identification.blocs), color = "black",size = 5) +
      ggplot2::ggtitle("Factorial representation of numeric variables") +
      ggplot2::xlab(paste(nom.comp[1], inertie[comp[1],2]," %")) +
      ggplot2::ylab(paste(nom.comp[2], inertie[comp[2],2]," %")) +
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
      ggplot2::theme_classic(base_size = 12) + ggplot2::annotate(
        geom = "segment",
        x = rep(0,nrow(weight.num)),
        xend = weight.num[,comp[1]],
        y = rep(0,nrow(weight.num)),
        yend = weight.num[,comp[2]],
        col = "black",
        arrow = ggplot2::arrow(length = grid::unit(0.2, "cm")),size = 0.70
      )  +
      ggplot2::annotate("path",
                        x = 0 + 1 * cos(seq(0, 2 * pi, length.out = 100)),
                        y = 0 + 1* sin(seq(0, 2 * pi, length.out = 100)),size = 0.80)

    #If there is a numeric supplementary variable
    if (supp.var == TRUE){
      loading.supp <- do.call(rbind.data.frame,res.MBPCAOS$coord.supp.num)

      facto.num <-
        facto.num  + ggplot2::annotate(
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
          col = "blue",size = 5
        )

    }

    return(facto.num)
  }

  #Blocs
  if (choice == "blocs"){
    cor.bloc <- data.frame(t(sapply(1:length(res.MBPCAOS$block.components),function(b)diag(cor(res.MBPCAOS$block.components[[b]],res.MBPCAOS$components)))))
    rownames(cor.bloc) <- names(res.MBPCAOS$block.components)

    graph.bloc <- ggplot2::ggplot(data = cor.bloc) +
      ggplot2::geom_label( ggplot2::aes(
        x = as.numeric(cor.bloc[,comp[1]]),
        y = as.numeric(cor.bloc[,comp[2]]),
        label = rownames(cor.bloc)
      ),
      size = 5) +
      ggplot2::ggtitle("Correlation of blocks components with global component") +
      ggplot2::xlab(paste(nom.comp[1], res.MBPCAOS$inertie[comp[1],2]," %")) +
      ggplot2::ylab(paste(nom.comp[2], res.MBPCAOS$inertie[comp[2],2]," %")) +
      ggplot2::geom_hline(
        yintercept = 0,
        linetype = "dotted",
        color = "black",
        size = 1
      ) + ggplot2::geom_vline(
        xintercept = 0,
        linetype = "dotted",
        color = "black",
        size = 1
      ) + ggplot2::theme_classic(base_size = 15)

    if(supp.var == TRUE){
      if(any(res.MBPCAOS$nature.supp == "nom" | res.MBPCAOS$nature.supp == "ord")){
        #Qualitative variables
        barycentre <-
          do.call(rbind.data.frame, res.MBPCAOS$coord.supp.quali)
        graph.bloc <-
          graph.bloc  +  ggplot2::annotate(
            geom = "label",
            x = barycentre[, comp[1]],
            y = barycentre[, comp[2]],
            label =  rownames(barycentre),
            size = 5,
            col = "blue"
          )
      }

      if(any(res.MBPCAOS$nature.supp == "num" )){
        #Numeric variables
        loading.supp <- do.call(rbind.data.frame,res.MBPCAOS$coord.supp.num)
        graph.bloc <-
          graph.bloc  + ggplot2::annotate(
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
            col = "blue",size = 5
          )
      }
    }


    return(graph.bloc)
  }

  #Squared loadings
  if (choice =="squared.loading"){
    sq.load  <- data.frame(t(data.frame(res.MBPCAOS$weights)))^2
    identification.blocs <- unlist(sapply(1:nb.bloc,function(x)rep(res.MBPCAOS$blocs.name[x],res.MBPCAOS$blocs[x],)))

    sq.load.graph <- ggplot2::ggplot(data = sq.load) +
      ggplot2::geom_label( ggplot2::aes(
        x = as.numeric(sq.load[,comp[1]]),
        y = as.numeric(sq.load[,comp[2]]),
        label = rownames(sq.load),
        fill = as.character(identification.blocs)
      ),size = 5) +
      ggplot2::xlab(paste(nom.comp[1], res.MBPCAOS$inertie[comp[1], 2], " %")) +
      ggplot2::ylab(paste(nom.comp[2], res.MBPCAOS$inertie[comp[2], 2], " %")) +
      ggplot2::geom_vline(
        xintercept = 0,
        linetype = "dotted",
        color = "black",
        size = 1
      ) +  ggplot2::geom_hline(
        yintercept = 0,
        linetype = "dotted",
        color = "black",
        size = 1
      )+
      ggplot2::theme_classic(base_size = 18) + ggplot2::guides(
        fill = ggplot2::guide_legend(
          title = "Blocs",
          override.aes =  ggplot2::aes(label = "")
        )
      ) + ggplot2::annotate(
        geom = "segment",
        x = rep(0,nrow(sq.load)),
        xend = sq.load[,comp[1]],
        y = rep(0,nrow(sq.load)),
        yend = sq.load[,comp[2]],
        col = "black",
        arrow = arrow(length = unit(0.2, "cm")),size = 0.70
      )
    return(sq.load.graph)
  }

  #Mixed
  if (choice == "mixed"){
    #NUM DATA
    data.graph.var <- data.frame(t(data.frame(res.MBPCAOS$weights)))
    weight.num <- data.graph.var[which(nature == "num"),]

    #NOM ET ORD DATA
    nb.modal <- unlist(lapply(category.coord,nrow))
    nb.var <- length(nb.modal)
    identification.variable <- as.vector(unlist(sapply(1:nb.var, function(j) {rep(colnames(res.MBPCAOS$data[,var.quali])[j],nb.modal[j])})))
    data.modal <- data.frame(data.modal,identification.variable)


    mixed.graph = ggplot2::ggplot() +
      ggplot2::geom_label(ggplot2::aes(
        x = as.numeric(weight.num[,comp[1]]),
        y = as.numeric(weight.num[,comp[2]]),
        label = rownames(weight.num)
      ),color = "black",size = 5) + ggplot2::annotate(
        geom = "segment",
        x = rep(0,nrow(weight.num)),
        xend = weight.num[,1],
        y = rep(0,nrow(weight.num)),
        yend = weight.num[,2],
        col = "black",
        arrow = arrow(length = unit(0.2, "cm")),size = 0.70
      ) +
      ggplot2::geom_hline(
        yintercept = 0,
        linetype = "dotted",
        color = "black",
        size = 1
      ) + ggplot2::geom_vline(
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
      size = 5) +
      ggplot2::geom_line(ggplot2::aes(
        x = as.numeric(data.modal[,2]),
        y = as.numeric(data.modal[,3]),
        group = as.character(identification.variable),
        col = as.character(identification.variable)
      ),
      size = 1,show.legend = F) +
      ggplot2::ggtitle("Mixed variables") +
      ggplot2::xlab(paste(nom.comp[1], res.MBPCAOS$inertie[comp[1],2]," %")) +
      ggplot2::ylab(paste(nom.comp[2], res.MBPCAOS$inertie[comp[2],2]," %")) +
      ggplot2::theme_classic(base_size = 18) + ggplot2::guides(fill = ggplot2::guide_legend(title = "Nominal variables",
                                                                 override.aes = ggplot2::aes(label = "")))
    return(mixed.graph)
  }

}
