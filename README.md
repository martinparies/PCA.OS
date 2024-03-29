# PCAOS

PCAOS package allows to perform Single and MultiBlock exploratory analysis of datasets composed of variables with different natures.

# How do you install an use the latest version of PCAOS available on GitHub?

```{r}
if (!require("devtools")) install.packages("devtools")
library(devtools)
install_github("martinparies/PCA.OS")
library(PCA.OS)
help(PCA.OS)
```

# 1. Single block data
```{r}
data("antibiotic")
help(antibiotic)
```

# 1.1 Setting nature of each variable
```{r}
#Manually
level.scale <- rep(NA,ncol(antibiotic)) #Setting level.scale argument
level.scale[c(3,4)] <- "num"
level.scale[c(6:14)] <- "nom"
level.scale[c(1,15)] <- "ord"

#Or using nature.variables()
level.scale <- rep(NA,ncol(antibiotic))
res.nature <- nature.variables(antibiotic)
level.scale [res.nature$p.numeric] <- "num"
level.scale [res.nature$p.quali] <- "nom"
#Warning; the ordinal nature of variables can not be detected automaticaly.
level.scale[c(1,15)] <- "ord"
```

# 1.2 Choice of number of component
```{r}
help(choice.component)
res.choice <- choice.component(antibiotic,level.scale)
res.choice
```

# 1.3 Single block exploratory analysis : PCAOS
```{r}
help(PCAOS)
res.PCAOS <- PCAOS(
  data = antibiotic,
  level.scale = level.scale,
  nb.comp = 2)
```

# 1.4 Plots
```{r}
help(plot.PCAOS)
#Individuals
plot.PCAOS(
  x = res.PCAOS,
  choice = "ind",
  coloring.indiv = antibiotic$Atb.conso,
  size.legend = 12,
  size.label = 4
)

# Variables
plot.PCAOS(x = res.PCAOS,choice = "all.var")
```

# 2. Multiblock data
```{r}
data('antibiotic')
antb.uses <- antibiotic[,c('Atb.conso','Atb.Sys')]
health <- antibiotic[,c('Age','Loss')]
vet.practices <- antibiotic[,c(6:15)]
antibiotic.MB <- data.frame(antb.uses,health,vet.practices)
```

# 2.1 Defining the blocks
```{r}
blocks.name =  c("antibiotic.uses","Health.of.turkeys","Veterinary.practices")
blocks <- c(2,2,10)
```

# 2.2 Level of scaling
```{r}
level.scale.MB <- rep(NA,ncol(antibiotic.MB))
res.nature <- nature.variables(antibiotic.MB)
level.scale.MB [res.nature$p.numeric] <- "num"
level.scale.MB [res.nature$p.quali] <- "nom"
#Warning; the ordinal nature of variables can not be detected automaticaly.
level.scale.MB[c(1,14)] <- "ord"
```

# 2.3 Choice of number of component
```{r}
help(choice.component.MB)
res.choice.MB <- choice.component.MB(antibiotic.MB,level.scale.MB, blocks , blocks.name, block.scaling = 'inertia')
res.choice.MB
```

# 2.4 Multiblock exploratory analysis : MBPCAOS
```{r}
res.MBPCAOS <- MBPCAOS(data = antibiotic.MB,
                       level.scale = level.scale.MB,
                       blocks = blocks,
                       blocks.name = blocks.name,
                       nb.comp = 3)

```

# 2.5 Blocks graph
```{r}
plot.MBPCAOS(res.MBPCAOS,choice = 'blocks')
```
