## Run GCM analysis, write model results

## Before: otoliths.csv, tags.csv (data)
## After:  gcm_results.RData (model)

library(TAF)
library(RTMB)
library(fishgrowth)

mkdir("model")

## Read data
otoliths <- read.taf("data/otoliths.csv")
tags <- read.taf("data/tags.csv")

## Estimate
fm <- lm(len~age, otoliths)
L0 <- coef(fm)[[1]]
rmax <- coef(fm)[[2]]

## Parameter list
par <- list(
  L0 = L0,
  log_rmax = log(rmax),
  log_k = log(1),
  t50 = 1,
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
map <- list(L0=factor(NA), log_rmax=factor(NA))  # fix L0 and rmax
model <- gcm(par, data, map=map)
fit <- nlminb(model$par, model$fn, model$gr,
              control=list(eval.max=1e4,iter.max=1e4))
report <- model$report()
sdreport <- sdreport(model, getReportCovariance=FALSE)

## Save results
save(model, fit, report, sdreport, file="model/gcm_results.RData")
