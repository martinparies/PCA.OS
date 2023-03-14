#' plot.MBPCAOS
#'
#' Visualisation of results from MBPCAOS method. See details for available plots.
#'
#' @param x an object of class MBPCAOS
#' @param choice the graph to plot possible values are "screeplot","quantif","indiv","cor","modalities","mixed","squared loadings". See Details.
#' @param comp a length 2 vector with the components to plot
#' @param supp.var boolean (FALSE by default), if TRUE supplementary variables are added in factorial representation
#' @param size.label size of label in graphs (all plots).
#' @param size.legend size of label in graphs (all plots).
#' @param sub.var.quantif a vector with variable of interest (quantification plots).
#' @param coloring.indiv a vector of length N to color individuals. If NULL, no coloring is applied (individuals plot).
#' @param ellipse boolean (FALSE by default), if TRUE, draw ellipses around categories of the qualitative variable considered as supplementary (individuals plot).
#' @param level.conf level of confidence ellipses (individuals plot).
#' @param min.contribution (all.var plot) Variables with a contribution (i.e loading) lower than this value will not be plotted in the 'all.var' graph (useful for dataset with a lot of variables)
#' @param label.cat if == 'var+cat', the name of the variable is included in the labels of the categories; if == 'cat', only the name of the categorie is plotted (name of categories should be unique) (qualitative and all.var plot)
#' @param ordinal.as.direction boolean (FALSE by default); if TRUE ordinal variables are represented as vectors, from the first categorie to the last one (qualitative and all.var plot).
#' @param label.size.freq boolean (FALSE by default); if TRUE size of categories are proportional to their citation frequencies (qualitative and all.var plot).
#' @param sub.bloc a scalar indicating the block to plot for variables graphs (i.e if sub.bloc == 1, variables of the first block are plotted) (all.var plot).
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
#'   \item  qualitative: factorial representation of qualitatives variables trough the representation of it's categories. Coordinates of each category is calculted such as the single quantification of the category multiplied by the loading of the associated variable.
#' }
#' For mixed  variables
#'  \itemize{
#'   \item  all.var: factorial representation of all variables (weight for numeric variables, and categories for qualitative variables)
#'   \item  squared.loadings: plot of the squared loadings of all the variables.
#' }
#' For blocs
#'  \itemize{
#'   \item  blocks: contribution of each block to principal components. The contribution of one block is calculated as the squared sum of the loading of the variables in the block, divided by the block scaling of the block.
#' }
#'
#' All graph are ggplot object
#'
#' @examples
#'
#' #'data('antibiotic')
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
#' @export plot.MBPCAOS
#' @export
#'
#'
plot.MBPCAOS <-function(
  x,
  choice,
  comp = c(1,2),
  coloring.indiv = NULL,
  supp.var = FALSE,
  sub.var.quantif = NULL,
  sub.bloc = NULL,
  size.label = 3.5,
  size.legend = 12,
  label.cat = 'var+cat',
  min.contribution = 0,
  ellipse = FALSE,
  level.conf = 0.95,
  ordinal.as.direction = FALSE,
  label.size.freq = FALSE,...) {

  res.MBPCAOS <- x
  if(!inherits(res.MBPCAOS,"MBPCAOS")) stop("Non convenient object")

  nb.comp <- ncol(res.MBPCAOS$Dimension.reduction$components)
  check.plot.MB.arg(choice,res.MBPCAOS$Quantification$level.scale,res.MBPCAOS$Supp.var$level.scale.supp,supp.var,comp,nb.comp)

  #Information about the results
  data <- res.MBPCAOS$Quantification$data
  level.scale <- res.MBPCAOS$Quantification$level.scale
  variables <- colnames(data)
  nb.var <- length(variables)

  #Blocks
  nb.block <- length(res.MBPCAOS$Blocks$blocks.name)
  blocks <- res.MBPCAOS$Blocks$blocks
  blocks.name <- res.MBPCAOS$Blocks$blocks.name

  #Dimension reduction
  components <- res.MBPCAOS$Dimension.reduction$components
  nb.comp <- ncol(components)
  weights <- res.MBPCAOS$Dimension.reduction$weights
  inertie <- res.MBPCAOS$Dimension.reduction$inertia

  #QUANTIFICATIONS
  quantification <- vector("list", length = nb.var)
  quantification.nom <- res.MBPCAOS$Quantification$quantification.categories.nom
  quantification.ord <- res.MBPCAOS$Quantification$quantification.categories.ord
  quantified.data <- res.MBPCAOS$Quantification$quantified.data
  var.nom <- names(quantification.nom)
  var.ord <- names(quantification.ord)
  var.num <- names(data[which(level.scale == "num")])

  limit = c(min(quantified.data), max(quantified.data))
  compteur <- 1
  for (i in which(level.scale == "nom")) {
    quantification[[i]] <- quantification.nom[[compteur]]
    compteur <- compteur + 1
  }
  compteur <- 1
  for (i in which(level.scale == "ord")) {
    quantification[[i]] <- quantification.ord[[compteur]]
    compteur <- compteur + 1
  }
  compteur <- 1
  for (i in which(level.scale == "num")){
    quantification[[i]] <- cbind(data[,i],quantified.data[,i])
    compteur <- compteur + 1
  }
  names (quantification) <- variables
  nom.comp <- c(paste("CP",comp[1],sep=""),paste("CP",comp[2],sep=""))

  #Var.supp
  var.supp <- res.MBPCAOS$Supp.var$var.supp
  level.scale.supp <- res.MBPCAOS$Supp.var$level.scale.supp
  coord.supp.num <- res.MBPCAOS$Supp.var$coord.supp.num
  coord.supp.quali <- res.MBPCAOS$Supp.var$coord.supp.quali

  #Screeplot
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

  #Quantification
  if (choice == "quantif"){
    graphs.list <- list(NA)
    modalities <- NULL
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
        colnames(quantification[[var]]) <- c("values","quantifications")
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

  #Individuals
  if (choice == "ind"){
    nom.indiv <- rownames(data)
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
      barycentre <-
        do.call(rbind.data.frame, res.MBPCAOS$Supp.var$barycenters)

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

  #Numeric variables
  if(choice == "all.var" & !any(level.scale == 'nom'| level.scale == 'ord')){choice <- 'numeric'}
  if(choice == "qualitative" & !any(level.scale == 'nom'| level.scale == 'ord')){choice <- 'numeric'}
  if (choice == "numeric"){
    data.graph.var <- data.frame(t(data.frame(weights)))
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

  #COMMON ELEMENTS
  if (choice == "qualitative" | choice == "all.var"){
    category.coord <- list(NULL)
    compteur <- 1
    var.quali <- which(level.scale == "nom" | level.scale =="ord")
    nb.var.quali <- length(var.quali)

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
          as.numeric(quantification[[var]]) * weights[[var]][comp[1]],
          as.numeric(quantification[[var]]) * weights[[var]][comp[2]]
        )
      colnames(category.coord[[compteur]]) <- c("Modalites",nom.comp)
      category.coord[[compteur]] <- data.frame(category.coord[[compteur]])
      compteur <- compteur + 1
    }
    names(category.coord) <- colnames(data[,var.quali])
    data.modal <- category.coord[[1]]
    if (length(category.coord)>1){
      for (i in 2:length(category.coord)){
        data.modal <- rbind(data.modal,category.coord[[i]])
      }
    }
    colnames(data.modal) <- colnames(category.coord[[1]])
    nb.modal <- unlist(lapply(category.coord,nrow))

    identification.blocs <- unlist(sapply(1:nb.block,function(x)rep(blocks.name[x],blocks[x],)))

    #Construction of matrix mixed, containing all coordinates of all variables and categories
    #quali var
    nb.modal <- unlist(lapply(category.coord,nrow))
    nb.var <- length(nb.modal)

    #Variables
    if(nb.var == 0){nb.modal = 0}
    i.variables <- as.vector(unlist(sapply(1:nb.var, function(j) {rep(colnames(data[,var.quali,drop = F])[j],nb.modal[j])})))

    #Blocs
    i.blocs <- unlist(sapply(1:nb.block,function(x)rep(blocks.name[x],blocks[x]),simplify = F))
    i.blocs.quali <- identification.blocs[which(level.scale == "nom" | level.scale =="ord")]
    i.blocs.num <- identification.blocs[which(level.scale == "num")]
    i.blocs.quali <- unlist(sapply(1:nb.var,function(x)rep(i.blocs.quali[x],nb.modal[x],),simplify = F))

    #Quali var
    nature.quali <- NULL
    level.scale.quali <- level.scale[which(level.scale == 'nom'|level.scale == 'ord')]
    for (i in 1:length(nb.modal)){
      nature.quali <- c(nature.quali,rep(level.scale.quali[i],nb.modal[i]))
    }
    p.nom <- which(nature.quali == 'nom')
    p.ord <- which(nature.quali == 'ord')
    data.modal <- cbind(data.modal,i.blocs.quali,i.variables,nature.quali)
    row.names(data.modal) <- data.modal[,1]

    #Numeric var
    if (any(level.scale == 'num')){
      p.num <- 1:length(which(level.scale == "num"))
      p.nom <- p.nom + length(p.num)
      weight.num <- data.frame(t(data.frame(weights)))
      weight.num <- weight.num[which(level.scale == "num"),]
      weight.num <- cbind(weight.num[,comp],i.blocs.num,'num')
      nb.var.num <- nrow(weight.num)
      colnames(weight.num) <- c(nom.comp[1],nom.comp[2],"i.blocs.num",'nature')
      #Mixed var
      nature <- c(weight.num[,4],data.modal[,6])
      na <- rep(NA,(sum(nb.modal) + nb.var))
      mixed <- data.frame()
      mixed <- rbind(weight.num[,c(1,2)],data.modal[,c(2,3)])
      mixed <- cbind(mixed,c(i.blocs.num,i.blocs.quali))
      mixed <- cbind(mixed,c(row.names(weight.num),i.variables),nature)
    }else{
      nature <- data.modal[,6]
      mixed <- data.frame()
      mixed <- rbind(data.modal[,c(2,3)])
      mixed <- cbind(mixed,i.blocs.quali)
      mixed <- cbind(mixed,i.variables,nature)
    }

    #Calcul frequence de citation
    data.quali <- data.frame(data[,var.quali,drop = F])
    freq <- unlist(sapply(1:ncol(data.quali),function(var){table(data.quali[,var])},simplify = F))
    freq <- freq/nrow(data) *100
    mixed[,6] <- c(rep(NA,length(which(mixed[,5] == "num"))),freq)
    data.modal <- data.frame(data.modal,freq)
  }

  #MODALITIES REPRESENTATION
  if(choice == "all.var" & !any(level.scale == 'num')){choice <- 'qualitative'}
  if (choice == "qualitative" ){
    nb.var <- length(nb.modal)
    graph.modalites <- ggplot2::ggplot(data = data.modal) +
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

    if(ordinal.as.direction == TRUE & any(level.scale == 'ord')){
      data.ord <- mixed[which(mixed[,5] == 'ord'),]
      ord.var <- unique(data.ord[,4])
      x.max <- x.min <- rep(NA,length(ord.var))
      y.max <- y.min <- rep(NA,length(ord.var))
      for (i in 1:length(ord.var)){
        data.var <- data.ord[which(data.ord[,4] == ord.var[i]),]

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

    if(any(level.scale.supp == "nom" ) | any(level.scale.supp == "ord" )){
      coord.supp.quali <- do.call(rbind.data.frame,res.MBPCAOS$Supp.var$coord.supp.quali)
      graph.modalites <-
        graph.modalites  + ggplot2::annotate(geom = "label",x = coord.supp.quali[,comp[1]], y = coord.supp.quali[,comp[2]],label =  rownames(coord.supp.quali),size = size.label,col = "blue")
    }


    return(graph.modalites)

  }

  #Blocs
  if (choice == "blocks"){
    #cor.bloc <- data.frame(t(sapply(1:length(res.MBPCAOS$block.components),function(b)diag(cor(res.MBPCAOS$block.components[[b]],res.MBPCAOS$components)))))

    contrib <- res.MBPCAOS$Blocks$contrib.blocks
    for (i in 1:ncol(contrib)){
      contrib[,i] <- contrib[,i] / sum(contrib[,i]) * 100

    }

    graph.bloc <- ggplot2::ggplot(data = contrib) +
      ggplot2::geom_label(ggplot2::aes(
        x = as.numeric(contrib[,comp[1]]),
        y = as.numeric(contrib[,comp[2]]),
        label = rownames(contrib)
      ), size = size.label) +
      ggplot2::ggtitle("blocks according to their contributions to components") +
      ggplot2::xlab(paste(nom.comp[1], inertie[comp[1],1]," %")) +
      ggplot2::ylab(paste(nom.comp[2], inertie[comp[2],1]," %")) +
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
    sq.load  <- data.frame(t(data.frame(weights)))^2
    identification.blocks <- unlist(sapply(1:nb.block,function(x)rep(blocks.name[x],blocks[x],)))

    sq.load.graph <- ggplot2::ggplot(data = sq.load) +
      ggplot2::geom_label( ggplot2::aes(
        x = as.numeric(sq.load[,comp[1]]),
        y = as.numeric(sq.load[,comp[2]]),
        label = rownames(sq.load),
        fill = as.character(identification.blocks)
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
          title = "Blocks",
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
    p.quali <- which(mixed[,5] == 'nom' | mixed[,5] == 'ord')
    title <- 'Factorial representation of all variables'
    if(!is.null(sub.bloc)){
      mixed <- mixed[which(mixed[,3] == blocks.name[sub.bloc]),]
      nature.bloc <- level.scale[which(mixed[,3] == blocks.name[sub.bloc])]
      p.num <- which(mixed[,5] == 'num')
      p.nom <- p.quali <- which(mixed$nature == 'nom' | mixed$nature == 'ord')
      title <- paste("Variables of block",blocks.name[sub.bloc])
    }

    if (min.contribution > 0){
      weights <- t(data.frame(weights))
      invisible.var <- intersect(which(abs(as.numeric(weights[,1])) < min.contribution), which(abs(as.numeric(weights[,2])) < min.contribution))
      invisible.var <- rownames(weights[invisible.var,])
      row.to.supp <- which(mixed[,4] %in% invisible.var)
      if(length(row.to.supp) > 0){
        mixed <- mixed[-row.to.supp,]
        p.num <- which(mixed$nature == 'num')
        p.quali <- which(mixed$nature == 'nom' | mixed$nature == 'ord')
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
        size = 1) +
      ggplot2::ggtitle(title) +
      ggplot2::xlab(paste(nom.comp[1], inertie[comp[1], 1], " %")) +
      ggplot2::ylab(paste(nom.comp[2], inertie[comp[2], 1], " %")) +
      ggplot2::theme_classic(base_size = size.legend) +
      ggplot2::guides(fill = ggplot2::guide_legend(title = "Block",override.aes = ggplot2::aes(label = "")))

    if(ordinal.as.direction == TRUE & any(level.scale == 'ord')){
      data.ord <- mixed[which(mixed[,5] == 'ord'),]
      ord.var <- unique(data.ord[,4])
      x.max <- x.min <- rep(NA,length(ord.var))
      y.max <- y.min <- rep(NA,length(ord.var))
      for (i in 1:length(ord.var)){
        data.var <- data.ord[which(data.ord[,4] == ord.var[i]),]

        x.max[i] <- as.numeric(data.var[nrow(data.var),1])
        x.min[i] <- as.numeric(data.var[1,1])
        y.max[i] <- as.numeric(data.var[nrow(data.var),2])
        y.min[i] <- as.numeric(data.var[1,2])
      }

      mixed.graph <- mixed.graph +
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
      mixed.graph <- mixed.graph +
        ggplot2::geom_label(ggplot2::aes(
          x = as.numeric(mixed[p.quali,1]),
          y =  as.numeric(mixed[p.quali,2]),
          label = rownames(mixed[p.quali,]),
          fill = mixed[p.quali,3],
          size = mixed[p.quali, 6]),
          color = "black")
    }else{
      mixed.graph <- mixed.graph +
        ggplot2::geom_label(ggplot2::aes(
          x = as.numeric(mixed[p.quali,1]),
          y =  as.numeric(mixed[p.quali,2]),
          label = rownames(mixed[p.quali,]),
          fill = mixed[p.quali,3]),
          size = size.label,
        )
    }

    if(supp.var == TRUE){
      if(any(level.scale.supp == "num" )){
        #Numeric variables
        loading.supp <- do.call(rbind.data.frame,coord.supp.num)
        mixed.graph <-
          mixed.graph  + ggplot2::annotate(
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
        coord.supp.quali <- do.call(rbind.data.frame,res.MBPCAOS$Supp.var$coord.supp.quali)
        mixed.graph <-
          mixed.graph  + ggplot2::annotate(geom = "label",x = coord.supp.quali[,comp[1]], y = coord.supp.quali[,comp[2]],label =  rownames(coord.supp.quali),size = size.label,col = "blue")
      }

    }
  }
  return(mixed.graph)

}

