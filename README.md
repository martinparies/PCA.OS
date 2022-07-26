# PCAOS

PCAOS package perform Principal Component Analysis with Optimal Scaling, used to deal with mixed data.

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
level.scale <- rep(NA,ncol(antibiotic)) #Setting level.scale argument
level.scale[c(2,3,4)] <- "num"
level.scale[c(1,5,6,7,8,9,10,11,12,13,14,15)] <- "nom"
level.scale[c(1,15)] <- "ord"
```

# Choice of number of components
```{r}
help(choice.component)
res.choice <- choice.component(antibiotic,nature)
res.choice
```

# PCAOS Analysis
```{r}
help(PCAOS)
res.PCAOS <-
  PCAOS(
    data = antibiotic,
    level.scale = level.scale,
    rank.restriction = "one",
    nb.comp = 4,
    supp.var = 1
  )
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

#Variables
plot.PCAOS(res.PCAOS = res.PCAOS,choice = "mixed")
```
