

#' load libraries
library(readxl)
library(data.table)

#' load dad
dat <- read_excel(path = "Detection_methods_Bilal_Manuscript_Metagen.xlsx")
setDT(dat)

#' check number of participants
dat$Patient_No  |> unique() |> length()


#' Make binary variables for 'detection'
dat[ , qPCR := qPCR_detection=="Detected" ]
dat[ , direct := Direct_sequencing=="Detected" ]
dat[ , culture := Culture_detection=="Detected" ]

#' Make Ca positive and negative datasets
datpositive <- dat[EPA_results=="Ca+"]
datnegative <- dat[EPA_results=="Ca-"]

#' Make a dataset with one row per sample, 
#' indicating whether the filtered or unfiltered sample was detected.
#' 
#' Possible this isn't correct, how do you identify the 'mother' sample?
#' 
datpositivesingle <- 
  datpositive[ , .(qPCR=max(qPCR), culture=max(culture), direct=max(direct)), 
               by=.(Patient_No)]


#' One way tabulations for each method

datpositivesingle[, table(direct)]
datpositivesingle[, table(culture)]
datpositivesingle[, table(qPCR)]

#' The number for detected by culture looks wrong compared to the paper (27 not 26), check it
#'
#' Which patients are not detected by culture
#'
datpositivesingle[culture==0, Patient_No]

#'
#' Which patients are detected by culture
#'
datpositivesingle[culture==1, Patient_No]

#'
#' These are the McNemar tests for each pair of methods
#' 
datpositivesingle[, table(direct, culture) |> mcnemar.test()]
datpositivesingle[, table(qPCR, culture) |> mcnemar.test()]
datpositivesingle[, table(qPCR, direct) |> mcnemar.test()]


# Redo table1
sumperc = function(x) sprintf("%d (%.2f)",sum(x),100*mean(x))

dat[  , .(number=.N,culture=sumperc(culture),
          qPCR=sumperc(qPCR),
          direct=sumperc(direct)) , by=.(EPA_results,Type,Bristol_scale)][
  order(EPA_results, Type, -Bristol_scale)] |> knitr::kable(align="r")

#' The 'culture' column doesn't quite agree with the paper for these elements
#'
#' I think there are four filtered samples detected by culture with bristol scale = 7, but table 1 says 5.
#' 
dat[EPA_results=="Ca+" & Type=="Filtered" & Bristol_scale==7]
#'
#' I think there are four filtered samples detected by culture with bristol scale = 5, but table 1 says 3.
#' 
dat[EPA_results=="Ca+" & Type=="Filtered" & Bristol_scale==5]


