## Run GCM analysis, write model results

## Before: otoliths.csv, tags.csv (data)
## After:  gcm_2_results.RData (model)

library(TAF)
library(RTMB)
library(fishgrowth)

mkdir("model")

## Read data
otoliths <- read.taf("data/otoliths.csv")
tags <- read.taf("data/tags.csv")

## Parameter list
par <- list(
  L0 = log(20),
  log_rmax = log(120),
  log_k = log(2),
  t50 = 0,
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
  Loto = otoliths$len
)

## Run model
model <- gcm(par, data)
fit <- nlminb(model$par, model$fn, model$gr,
              control=list(eval.max=1e4,iter.max=1e4))
report <- model$report()
sdreport <- sdreport(model, getReportCovariance=FALSE)

## Save results
save(model, fit, report, sdreport, file="model/gcm_2_results.RData")
