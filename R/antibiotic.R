#' @title Antibiotic data
#'
#' @description This data set contains information about 15 variables about medical practice, antibiotic consumption
#' in 122 agricultural explotation of Turkey
#'
#' @format A data frame with 122 row (agricultural explotation) and 15 columns (variables).
#' \itemize{
#'   \item Atb.conso (Ordinal): Antibiotic consumption (low / medium / high)
#'   \item Duration (Numeric): Rearing time (in days)
#'   \item Age (Numeric): Average age at slaughter of turkeys (in days)
#'   \item Loss (Numeric): Percentage of turkeys that died during rearing
#'   \item Atb.Sys (Nominal): Systematic use of antibiotics (yes/no)
#'   \item Chick.an (Nominal): Systematic analysis of chicks (salmo, Coli) (yes/rarely-no)
#'   \item Autops (Nominal): Autopsy during rearing (yes/no)
#'   \item Vac (Nominal): change of vaccination plan during the year (yes/no)
#'   \item Med.coop (Nominal): Purchase of drugs from the cooperative(yes/no)
#'   \item Who.atb (Nominal): Who decides on antibiotic treatment? (farmer/veterinarian)
#'   \item Proced.pb (Nominal): What to do in case of problems (farmer/technician/lab/other)
#'   \item Duration.san (Nominal): Duration of sanitary vacuum (<3 weeks/> 3 weeks)
#'   \item Delivery.atb (Nominal): System for administering antibiotics in drinking water (tank, pump)
#'   \item Alim.atb (Nominal): Medicated feeding for curative purposes (yes/no)
#'   \item Agri.area (Ordinal): Useful agricultural area of the farm (<30ha / 30=<<50ha / >=50ha )
#' }
#'
#' @example
#' #Charging the data
#' data("antibiotic")
#' summary(antibiotic)
#' #Construction of the "nature" argument for this dataset
#' nature <- rep(NA,ncol(antibiotic)) #Setting nature argument
#' nature[c(2,3,4)] <- "num"
#' nature[c(1,5,6,7,8,9,10,11,12,13,14,15)] <- "nom"
#' nature[c(1,15)] <- "ord"
#' nature
#'
"antibiotic"
