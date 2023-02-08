# PCAOS

PCAOS package allows to perform Single and MultiBlock exploratory analysis of datasets composed of variables with different natures.
How do you install an use the latest version of PCAOS available on GitHub?


```{r}
if (!require("devtools")) install.packages("devtools")
library(devtools)
install_github("martinparies/PCA.OS")
library(PCA.OS)
help(PCA.OS)
```

# Data
```{r}
data("antibiotic")
help(antibiotic)
```

# Setting nature of each variable
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

# Choice of number of components
```{r}
help(choice.component)
res.choice <- choice.componentlevel.scale
res.choice
```

# Single block exploratory analysis : PCAOS
```{r}
help(PCAOS)
res.PCAOS <- PCAOS(
    data = antibiotic,
    level.scale = level.scale,
    rank.restriction = "one",
    nb.comp = 2)
```

# Plots
```{r}
help(plot.PCAOS)
#Individuals
plot.PCAOS(
  res.PCAOS = res.PCAOS,
  choice = "ind",
  coloring.indiv = antibiotic$Atb.conso,
  supp.var = TRUE,
  ellipse = TRUE,
  size.legend = 12,
  size.label = 4
)

# Variables
plot.PCAOS(res.PCAOS = res.PCAOS,choice = "all.var")
```

# Multiblock data
```{r}
data('antibiotic')
antb.uses <- antibiotic[,c('Atb.conso','Atb.Sys')]
health <- antibiotic[,c('Age','Loss')]
vet.practices <- antibiotic[,c(6:15)]
antibiotic.MB <- data.frame(antb.uses,health,vet.practices)
 ```

# Defining the blocks
 ```{r}
blocs.name =  c("antibiotic.uses","Health.of.turkeys","Veterinary.practices")
blocs <- c(2,2,10)
```

# Level of scaling
```{r}
level.scale.MB <- rep(NA,ncol(antibiotic.MB))
res.nature <- nature.variables(antibiotic.MB)
level.scale.MB [res.nature$p.numeric] <- "num"
level.scale.MB [res.nature$p.quali] <- "nom"
#Warning; the ordinal nature of variables can not be detected automaticaly.
level.scale.MB[c(1,14)] <- "ord"
```

# Choice of number of components
```{r}
help(choice.component.MB)
res.choice.MB <- choice.component.MB(antibiotic.MB,level.scale.MB, blocs , blocs.name, block.scaling = 'inertia')
res.choice.MB
```

Multiblock exploratory analysis : MBPCAOS
```{r}
res.MBPCAOS <- MBPCAOS(data = antibiotic.MB,
                       level.scale = level.scale.MB,
                       blocs = blocs,
                       blocs.name = blocs.name,
                       nb.comp = 3)

```

# Blocks graphs
```{r}
plot.MBPCAOS(res.MBPCAOS,choice = 'blocs')
```
