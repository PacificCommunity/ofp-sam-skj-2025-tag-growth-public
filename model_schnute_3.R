## Run Schnute case 3 analysis, write model results

## Before: otoliths.csv, tags.csv (data)
## After:  schnute_3_results.RData (model)

library(TAF)
library(RTMB)
library(fishgrowth)

mkdir("model")

## Read data
otoliths <- read.taf("data/otoliths.csv")
tags <- read.taf("data/tags.csv")

## Parameter list
par <- list(
  log_L1 = log(30),
  log_L2 = log(80),
  b = 1.0,
  log_sigma_min = log(1),
  log_sigma_max = log(2),
  log_age = log(tags$lenRel/60)
)

## Data list
data <- list(
  Lrel = tags$lenRel,
  Lrec = tags$lenRec,
  liberty = tags$libertyYears,
  Aoto = otoliths$age,
  Loto = otoliths$len,
  t1 = 0,
  t2 = 4
)

## Run model
model <- schnute3(par, data)
fit <- nlminb(model$par, model$fn, model$gr,
              control=list(eval.max=1e4,iter.max=1e4))
report <- model$report()
sdreport <- sdreport(model, getReportCovariance=FALSE)

## Save results
save(model, fit, report, sdreport, file="model/schnute_3_results.RData")
