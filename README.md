# PCAOS

PCAOS package perform Principal Component Analysis with Optimal Scaling, used to deal with mixed data.

How do you install an use the latest version of PCAOS available on GitHub?

```{r}
if (!require("devtools")) install.packages("devtools")
library(devtools)
install_github("martinparies/PCAOS")
library(PCA.OS)
help(PCA.OS)
```

# Data
```{r}
data("antibiotic")
help(antibiotic)
```

# Setting of nature of variables
```{r}
nature <- rep(NA,ncol(antibiotic))
help(PCA.OS::mix.data)
res.mix <- PCA.OS::mix.data(antibiotic)
nature[res.mix$p.numeric] <- "num"
nature[res.mix$p.quali] <- "nom"
#No automatic solution for ordinal data, setting by hand
nature[c(1,15)] <- "ord"
```

# Choice of components
```{r}
help(choice.component)
res.choice <- PCAOS::choice.component(antibiotic,nature)
res.choice
```

# Analysis and plot
```{r}
help(PCAOS)
res.PCAOS <-
  PCAOS::PCAOS(
    data = antibiotic,
    nature = nature,
    rank.restriction = "one",
    nb.comp = 4,
    supp.var = 1
  )

#PLOTS
help(plot.PCAOS)
#Individuals
PCAOS::plot.PCAOS(
  res.PCAOS = res.PCAOS,
  choice = "ind",
  coloring.indiv = antibiotic$Atb.conso,
  supp.var = TRUE,
  conf.ellipsises = TRUE,
  size.legend = 12,
  size.label = 4
)
#Variables
PCAOS::plot.PCAOS(res.PCAOS = res.PCAOS,choice = "mixed")
```
