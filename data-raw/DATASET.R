## code to prepare `DATASET` dataset goes here
#One block of mixed data
setwd("C:/Users/martin.paries/Desktop/These/Sript R/Package/PCAOSproject/data-raw")
antibiotic <- data.frame(readRDS("DataSel.rds"))
usethis::use_data(antibiotic, overwrite = TRUE)

#MB of mixed data
library(PCAmixdata)
data("gironde")
usethis::use_data(gironde, overwrite = TRUE)
