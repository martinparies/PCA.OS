#' plot.MBPCAOS
#'
#' Visualisation of results from MBPCAOS method. See details for available plot.
#'
#' @param res.MBPCAOS an object of class MBPCAOS
#' @param choice the graph to plot possible values are "screeplot","quantif","indiv","cor","modalities","mixed","squared loadings". See Details.
#' @param comp a length 2 vector with the components to plot
#' @param coloring.indiv A vector of length N to color individuals. If NULL, no coloring is applied.
#' @param sub.var.quantif a vector with variable of interest
#' @param supp.var boolean (FALSE by default), if TRUE supplementary variables are added in factorial representation
#' @param min.contribution variables with a contribution (i.e loading) lower than this value will not be plotted in the 'all.var' graph (useful for dataset with a lot of variables)
#' @param label.cat If == 'var+cat', the name of the variable is included in the labels of the categories; if == 'cat', only the name of the categorie is plotted (name of categories should be unique)
#' @param sub.bloc A scalar indicating the block to plot for variables graphs, works only with 'all.var' graphs (i.e if sub.bloc == 1, variables of the first block are plotted)
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
#'   \item  all.var: factorial representation of all variables (weight for numeric variables, and categories for qualitative variables)
#'   \item  squared.loadings: plot of the squared loadings of all the variables.
#' }
#' For blocs
#'  \itemize{
#'   \item  blocs: contribution of each block to principal components. The contribution of one block is calculated as the squared sum of the loading of the variables in the block, divided by the block scaling of the block.
#' }
#'
#' All graph are ggplot object
#'
#' @export plot.MBPCAOS
#' @export
#'
#'
plot.MBPCAOS <-function(
  res.MBPCAOS,
  choice,
  comp = c(1,2),
  coloring.indiv = NULL,
  supp.var = FALSE,
  sub.var.quantif = NULL,
  sub.bloc = NULL,
  size.label = 3.5,
  size.legend = 12,
  label.cat = 'var+cat',
  min.contribution = 0) {

  nb.comp <- ncol(res.MBPCAOS$components)
  check.plot.MB.arg(choice,res.MBPCAOS$level.scale,res.MBPCAOS$level.scale.supp,supp.var,comp,nb.comp)

  #Information about the results
  data <- res.MBPCAOS$data
  level.scale <- res.MBPCAOS$level.scale
  variables <- colnames(res.MBPCAOS$data)
  nb.comp <- ncol(res.MBPCAOS$T.tot)
  nb.var <- length(variables)
  nb.bloc <- length(res.MBPCAOS$blocs.name)
  quantification <- vector("list", length = nb.var)
  quantification.nom <- res.MBPCAOS$quantification.categories.nom
  quantification.ord <- res.MBPCAOS$quantification.categories.ord
  var.nom <- names(quantification.nom)
  var.ord <- names(quantification.ord)
  var.num <- names(res.MBPCAOS$data[which(level.scale == "num")])
  inertie <- res.MBPCAOS$inertia
  limit = c(min(res.MBPCAOS$quantified.data), max(res.MBPCAOS$quantified.data))

  compteur <- 1
  for (i in which(level.scale == "nom")) {
    quantification[[i]] <- res.MBPCAOS$quantification.categories.nom[[compteur]]
    compteur <- compteur + 1
  }
  compteur <- 1
  for (i in which(level.scale == "ord")) {
    quantification[[i]] <- res.MBPCAOS$quantification.categories.ord[[compteur]]
    compteur <- compteur + 1
  }
  compteur <- 1
  for (i in which(level.scale == "num")){
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
      ggplot2::theme_classic(base_size = size.legend)

    return (screeplot)
  }

  #Quantification
  if (choice == "quantif"){
    graphs.list <- list(NA)
    limit = c(min(res.MBPCAOS$quantified.data)-0.1, max(res.MBPCAOS$quantified.data)+0.1)
    for (var in 1:nb.var){
      if(level.scale[var] == "nom"){
        data.pour.graph <- cbind(as.numeric(quantification[[var]]),rownames(quantification[[var]]))
        colnames(data.pour.graph) <- c("quantification","modalities")
        graphs.list[[var]] <- ggplot2::ggplot(data = data.frame(data.pour.graph),
                                              ggplot2::aes(x=modalities, y=as.numeric(quantification))) +
          ggplot2::geom_point(size = 3) +
          #ggplot2::scale_y_continuous(limits=limit) +
          ggplot2::scale_x_discrete(limits = levels(data[,var])) +
          ggplot2::xlab("Categories") +
          ggplot2::ylab("Quantifications") +
          ggplot2::ggtitle(variables[var])+
          ggplot2::theme_classic(base_size = size.legend)
      }
      if (level.scale[var] == "ord"){
        data.pour.graph <- cbind(as.numeric(quantification[[var]]),rownames(quantification[[var]]))
        order.modal <- rownames(quantification[[var]])
        #data.pour.graph <- data.frame(quantification[[var]])
        colnames(data.pour.graph) <- c("quantification","modalities")
        graphs.list[[var]] <- ggplot2::ggplot(data = data.frame(data.pour.graph), ggplot2::aes(x=modalities, y=as.numeric(quantification))) +
          ggplot2::geom_point(size = 3) +
          #ggplot2::scale_y_continuous(limits=limit) +
          ggplot2::scale_x_discrete(limits = levels(data[,var])) +
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

  #Individuals
  if (choice == "ind"){
    nom.indiv <- rownames(res.MBPCAOS$data)
    limit.x = c(min(res.MBPCAOS$components[,comp[1]]), max(res.MBPCAOS$components[,comp[1]]))
    limit.y = c(min(res.MBPCAOS$components[,comp[2]]), max(res.MBPCAOS$components[,comp[2]]))

    Graph.observations <-
      ggplot2::ggplot(
        data.frame(res.MBPCAOS$components),
        ggplot2::aes(
          x = res.MBPCAOS$components[,comp[1]],
          y = res.MBPCAOS$components[,comp[2]],
          label = rownames(res.MBPCAOS$components)
        )
      ) +
      ggplot2::ggtitle("Factorial representation of individuals") +
      ggplot2::geom_point() +
      ggplot2::geom_text(ggplot2::aes(label=rownames(res.MBPCAOS$components), color = coloring.indiv),hjust=1, vjust=1,size = size.label)+
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

    if(choice == "qualitative" & supp.var == TRUE) {message("Modalities are represented as barycenter")}

    #Si il y a une variable qualitative supplementaire
    if (supp.var == TRUE){
      barycentre <-
        do.call(rbind.data.frame, res.MBPCAOS$coord.supp.quali)

      Graph.observations <-
        Graph.observations  + ggplot2::annotate(geom = "label",x = barycentre[,comp[1]], y = barycentre[,comp[2]],label =  rownames(barycentre),size = size.label,col = "blue")

      if(ellipse == TRUE){
        mat = list(NULL)

        for (var.supp in 1:ncol(res.MBPCAOS$quali.var.supp)){

          var.quali = res.MBPCAOS$quali.var.supp[,var.supp]
          coord.ellipse <- cbind.data.frame(ellipses.coord(var.quali = var.quali,components = res.MBPCAOS$components,level.conf = level.conf))
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

    return(Graph.observations)  }

  #COMMON ELEMENTS
  if (choice == "qualitative" | choice == "all.var"){
    category.coord <- list(NULL)
    compteur <- 1
    var.quali <- which(level.scale == "nom" | level.scale =="ord")
    nb.var.quali <- length(var.quali)

    #if (res.MBPCAOS$summary$rank == "one"){
    #Computation of categories coordinates
    if (choice == "qualitative" | choice == "all.var"){
      weigths <- res.MBPCAOS$weights
    }

    for (var in var.quali){
      category.coord[[compteur]] <- rep(NA,3)
      if(label.cat == 'var+cat'){
        labels.of.cat <- paste(colnames(data[,var,drop = FALSE]),rownames(quantification[[var]]),sep = "_")
      }
      if(label.cat == 'cat'){
        labels.of.cat <- paste(rownames(quantification[[var]]),sep = "_")
      }
      category.coord[[compteur]] <- rep(NA,3)
      category.coord[[compteur]] <-
        cbind(
          labels.of.cat,
          #paste(rownames(quantification[[var]]),sep = "_"),
          as.numeric(quantification[[var]]) * res.MBPCAOS$weights[[var]][comp[1]],
          as.numeric(quantification[[var]]) * res.MBPCAOS$weights[[var]][comp[2]]
        )
      colnames(category.coord[[compteur]]) <- c("Modalites",nom.comp)
      category.coord[[compteur]] <- data.frame(category.coord[[compteur]])
      compteur <- compteur + 1
    }
    names(category.coord) <- colnames(res.MBPCAOS$data[,var.quali])
    data.modal <- category.coord[[1]]
    if (length(category.coord)>1){
      for (i in 2:length(category.coord)){
        data.modal <- rbind(data.modal,category.coord[[i]])
      }
    }
    colnames(data.modal) <- colnames(category.coord[[1]])
    nb.modal <- unlist(lapply(category.coord,nrow))

    #Identification of the number of modalities per variable (for fill argument)
    # if (coloring.cat == "variables"){
    #   fill.arg <- as.vector(unlist(sapply(1:nb.var, function(j) {rep(colnames(res.MBPCAOS$data[,var.quali])[j],nb.modal[j])})))
    #   legend <- "Variables"
    # }
    # if (coloring.cat == "blocs"){
    #   identification.blocs <- unlist(sapply(1:nb.bloc,function(x)rep(res.MBPCAOS$blocs.name[x],res.MBPCAOS$blocs[x],)))
    #   identification.blocs <- identification.blocs[which(nature == "nom" | nature =="ord")]
    #   fill.arg <- unlist(sapply(1:nb.var.quali,function(x) {rep(identification.blocs[x],nb.modal[x])},simplify = F))
    #   legend <- "Blocs"
    # }

    identification.blocs <- unlist(sapply(1:nb.bloc,function(x)rep(res.MBPCAOS$blocs.name[x],res.MBPCAOS$blocs[x],)))

    ## Detection of limits of graphs
    # max.x <- as.numeric(max(data.modal[,2]))
    # min.x <- as.numeric(min(data.modal[,2]))
    #
    # max.y <- as.numeric(max(data.modal[,3]))
    # min.y <- as.numeric(min(data.modal[,3]))
    #
    # limit.x <- c(min.x-0.2,max.x+0.2)
    # limit.y <- c(min.y-0.2,max.y+0.2)

    #Construction of matrix mixed, containing all coordinates of all variables and categories
    #quali var
    nb.modal <- unlist(lapply(category.coord,nrow))
    nb.var <- length(nb.modal)

    #Variables
    if(nb.var == 0){nb.modal = 0}
    i.variables <- as.vector(unlist(sapply(1:nb.var, function(j) {rep(colnames(res.MBPCAOS$data[,var.quali,drop = F])[j],nb.modal[j])})))
    #data.modal <- cbind(data.modal,fill.arg)

    #Blocs
    i.blocs <- unlist(sapply(1:nb.bloc,function(x)rep(res.MBPCAOS$blocs.name[x],res.MBPCAOS$blocs[x],),simplify = F))
    i.blocs.quali <- identification.blocs[which(level.scale == "nom" | level.scale =="ord")]
    i.blocs.num <- identification.blocs[which(level.scale == "num")]
    i.blocs.quali <- unlist(sapply(1:nb.var,function(x)rep(i.blocs.quali[x],nb.modal[x],),simplify = F))

    #Quali var
    p.nom <- 1:sum(nb.modal)
    data.modal <- cbind(data.modal,i.blocs.quali,i.variables,'quali')
    row.names(data.modal) <- data.modal[,1]

    #Numeric var
    p.num <- 1:length(which(level.scale == "num"))
    p.nom <- p.nom + length(p.num)
    if(choice == "qualitative" | choice == "all.var"){
      weight.num <- data.frame(t(data.frame(res.MBPCAOS$weights)))
      weight.num <- weight.num[which(level.scale == "num"),]
      weight.num <- cbind(weight.num[,comp],i.blocs.num,'num')
    }
    nb.var.num <- nrow(weight.num)
    colnames(weight.num) <- c(nom.comp[1],nom.comp[2],"i.blocs.num",'nature')

    #Mixed var
    nature <- c(weight.num[,4],data.modal[,6])
    na <- rep(NA,(sum(nb.modal) + nb.var))
    mixed <- data.frame()
    mixed <- rbind(weight.num[,c(1,2)],data.modal[,c(2,3)])
    mixed <- cbind(mixed,c(i.blocs.num,i.blocs.quali))
    mixed <- cbind(mixed,c(row.names(weight.num),i.variables),nature)

    #Calcul frequence de citation
    data.quali <- data.frame(res.MBPCAOS$data[,var.quali,drop = F])
    freq <- unlist(sapply(1:ncol(data.quali),function(var){table(data.quali[,var])},simplify = F))
    freq <- freq/nrow(res.MBPCAOS$data) *100

    mixed[,6] <- c(rep(NA,length(which(mixed[,5] == "num"))),freq)

  }

  #MODALITIES REPRESENTATION
  if (choice == "qualitative" ){
    nb.var <- length(nb.modal)

    graph.modalites <- ggplot2::ggplot(data = data.modal) +
      ggplot2::geom_label(ggplot2::aes(
        x = as.numeric(data.modal[,2]),
        y =  as.numeric(data.modal[,3]),
        label = data.modal[, 1],
        size = freq
      )) +
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
    # graph.modalites <- graph.modalites + ggplot2::geom_line(ggplot2::aes(
    #   x = as.numeric(data.modal[,2]),
    #   y = as.numeric(data.modal[,3]),
    #   group = i.variables,
    #   col = as.character(i.variables)
    # ),
    # size = 1,show.legend = F)

    return(graph.modalites)

  }

  #Numeric variables
  if(choice == 'all.var' & !any(level.scale == 'nom' |level.scale=='ord')){choice = 'numeric'}
  if (choice == "numeric"){
    weight <- data.frame(t(data.frame(res.MBPCAOS$weights)))
    weight.num <- weight[which(level.scale == "num"),]
    identification.blocs <- unlist(sapply(1:nb.bloc,function(x)rep(res.MBPCAOS$blocs.name[x],res.MBPCAOS$blocs[x],)))
    identification.blocs <- identification.blocs[which(level.scale == "num")]

    graph.num <-
      ggplot2::ggplot(data = data.frame(weight.num),
                      ggplot2::aes (
                        x = as.numeric(weight.num[,comp[1]]),
                        y = as.numeric(weight.num[,comp[2]])
                      )) +
      ggplot2::geom_label(ggplot2::aes(label = row.names(weight.num),fill = identification.blocs), color = "black",size = size.label) +
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
        arrow = ggplot2::arrow(length = grid::unit(0.2, "cm")),size = 0.70
      )  +
      ggplot2::annotate("path",
                        x = 0 + 1 * cos(seq(0, 2 * pi, length.out = 100)),
                        y = 0 + 1* sin(seq(0, 2 * pi, length.out = 100)),size = 0.80) +
      ggplot2::guides(fill = ggplot2::guide_legend(title = "Blocs", override.aes = ggplot2::aes(label = "")))

    #If there is a numeric supplementary variable
    if (supp.var == TRUE){
      loading.supp <- do.call(rbind.data.frame,res.MBPCAOS$coord.supp.num)

      graph.num <-
        graph.num  + ggplot2::annotate(
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

    return(graph.num)
  }

  #Blocs
  if (choice == "blocks"){
    # cor.bloc <- data.frame(t(sapply(1:length(res.MBPCAOS$block.components),function(b)diag(cor(res.MBPCAOS$block.components[[b]],res.MBPCAOS$components)))))
    nb.block <- length(res.MBPCAOS$blocs)
    blocs <- res.MBPCAOS$blocs
    index.group <- unlist(sapply(1:nb.block,function(i) rep(i,res.MBPCAOS$blocs[i])))
    blocs.list <- list(NULL)
    cumul <- c(1,sapply(1:nb.block,function(j){sum(blocs[1:j])}))
    cumul2 <- c(1,sapply(2:length(cumul),function(x)cumul[x]+1))
    blocs.list <- sapply(1:nb.block,function(j){cumul2[j]:cumul[j+1]},simplify = FALSE)
    weight <- t(data.frame(res.MBPCAOS$weight))
    for (i in 1:ncol(weight)){
      weight[,i] <- (weight[,i]^2) / sum(weight[,i]^2)*100
    }
    weight.list <- sapply(1:nb.block,function(j) (weight[blocs.list[[j]],]) / (res.MBPCAOS$block.scaling[j]) )

    contrib <- matrix(NA,nb.block,ncol(res.MBPCAOS$components))
    for (i in which(unlist(lapply(weight.list,is.matrix)))) {
      weight.list[[i]] <- apply(weight.list[[i]],2,sum)
    }

    contrib <- data.frame(t(data.frame(weight.list)))
    rownames(contrib) <- blocs.name
    rownames(contrib) <- names(res.MBPCAOS$block.components)

    graph.bloc <- ggplot2::ggplot(data = contrib) +
      ggplot2::geom_label(ggplot2::aes(
        x = as.numeric(contrib[,comp[1]]),
        y = as.numeric(contrib[,comp[2]]),
        label = rownames(contrib)
      ), size = size.label) +
      ggplot2::ggtitle("Blocks according to their contributions to components") +
      ggplot2::xlab(paste(nom.comp[1], res.MBPCAOS$inertia[comp[1],1]," %")) +
      ggplot2::ylab(paste(nom.comp[2], res.MBPCAOS$inertia[comp[2],1]," %")) +
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
      ),size = size.label) +
      ggplot2::xlab(paste(nom.comp[1], res.MBPCAOS$inertie[comp[1], 1], " %")) +
      ggplot2::ylab(paste(nom.comp[2], res.MBPCAOS$inertie[comp[2], 1], " %")) +
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
      ggplot2::theme_classic(base_size = size.legend) + ggplot2::guides(
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
        arrow = ggplot2::arrow(length = grid::unit(0.2, "cm")),size = 0.70
      )
    return(sq.load.graph)
  }

  #Mixed
  if (choice == "all.var"){
    blocs.name <- names(res.MBPCAOS$block.components)
    if(!is.null(sub.bloc)){
      mixed <- mixed[which(mixed[,3] == blocs.name[sub.bloc]),]
      nature.bloc <- level.scale[which(mixed[,3] == blocs.name[sub.bloc])]
      p.num <- which(mixed[,5] == 'num')
      p.nom <- which(mixed[,5] == 'quali')

      mixed.graph <- ggplot2::ggplot() +
        ggplot2::geom_label(ggplot2::aes(
          x = as.numeric(mixed[p.num,1]),
          y = as.numeric(mixed[p.num,2]),
          label = rownames(mixed[p.num,]),
        ),size = size.label) + ggplot2::annotate(
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
        ) + ggplot2::geom_vline(
          xintercept = 0,
          linetype = "dotted",
          color = "black",
          size = 1
        ) +
        ggplot2::geom_label(ggplot2::aes(
          x = as.numeric(mixed[p.nom,1]),
          y =  as.numeric(mixed[p.nom,2]),
          label = row.names(mixed[p.nom,]),
          size = mixed[p.nom,6]
        )) +
        ggplot2::ggtitle(paste("Variables of blocs",blocs.name[sub.bloc])) +
        ggplot2::xlab(paste(nom.comp[1], res.MBPCAOS$inertia[comp[1], 1], " %")) +
        ggplot2::ylab(paste(nom.comp[2], res.MBPCAOS$inertia[comp[2], 1], " %")) +
        ggplot2::theme_classic(base_size = size.legend) +
        ggplot2::guides(fill = ggplot2::guide_legend(title = "Blocs",override.aes = ggplot2::aes(label = "")))

    }else{

      # data.ord <- mixed[which(mixed[,4] == 'ord'),]
      # ord.var <- unique(data.ord[,3])
      # x.max <- x.min <- rep(NA,length(ord.var))
      # y.max <- y.min <- rep(NA,length(ord.var))
      # for (i in 1:length(ord.var)){
      #   data.var <- data.ord[which(data.ord[,3] == ord.var[i]),]
      #
      #   x.max[i] <- as.numeric(data.var[nrow(data.var),1])
      #   x.min[i] <- as.numeric(data.var[1,1])
      #   y.max[i] <- as.numeric(data.var[nrow(data.var),2])
      #   y.min[i] <- as.numeric(data.var[1,2])
      # }

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

      mixed.graph <- ggplot2::ggplot() +
        ggplot2::geom_label(ggplot2::aes(
          x = as.numeric(mixed[p.num,1]),
          y = as.numeric(mixed[p.num,2]),
          label = rownames(mixed[p.num,]),
          fill = mixed[p.num,3]
        ),size = size.label) + ggplot2::annotate(
          geom = "segment",
          x = rep(0,nrow(mixed[p.num,])),
          xend = as.numeric(mixed[p.num,1]),
          y = rep(0,nrow(mixed[p.num,])),
          yend = as.numeric(mixed[p.num,2]),
          col = "black",
          arrow = ggplot2::arrow(length = grid::unit(0.2, "cm")),size = 0.50
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
        # ) + ggplot2::annotate(
        # geom = "segment",
        # x = x.min,
        # xend = x.max,
        # y = y.min,
        # yend = y.max,
        # col = "black",
        # arrow = ggplot2::arrow(length = grid::unit(0.45, "cm")),size = 0.70
      ) + ggplot2::geom_label(ggplot2::aes(
          x = as.numeric(mixed[p.nom,1]),
          y =  as.numeric(mixed[p.nom,2]),
          label = row.names(mixed[p.nom,]),
          fill = as.character(mixed[p.nom,3])
          #,size = mixed[p.nom,6]
        ),size = size.label) +
        # ggplot2::geom_line(ggplot2::aes(
        #   x = as.numeric(mixed[p.nom,1]),
        #   y = as.numeric(mixed[p.nom,2]),
        #   group = mixed[p.nom,4]
        #   #(col = as.character(mixed[p.nom,3])
        #   ),
        # size = 0.5,show.legend = F,color = "black") +
        ggplot2::ggtitle("All variables") +
        ggplot2::xlab(paste(nom.comp[1], res.MBPCAOS$inertia[comp[1], 1], " %")) +
        ggplot2::ylab(paste(nom.comp[2], res.MBPCAOS$inertia[comp[2], 1], " %")) +
        ggplot2::theme_classic(base_size = size.legend) +
        ggplot2::guides(fill = ggplot2::guide_legend(title = "Blocs",override.aes = ggplot2::aes(label = "")))

    }


    return(mixed.graph)

  }

}
