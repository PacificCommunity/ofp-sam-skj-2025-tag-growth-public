## Run GCM analysis, write model results

## Before: gcm.cpp (boot/software),
##         otoliths.csv, tags.csv (data)
## After:  gcm_1_results.RData (model)

library(TAF)
library(TMB)

mkdir("model")

## Read data
otoliths <- read.taf("data/otoliths.csv")
tags <- read.taf("data/tags.csv")

## Specify Lfix
Lfix <-  34  # could use 58

## Parameter list
par <- list(
  log_Linf = log(73),
  log_rmax = log(20),
  log_k = log(0.6),
  log_Lfix = log(Lfix),
  log_sigma_1 = log(0.0000001),
  log_sigma_2 = log(2.3),
  log_age = log(tags$lenRel/Lfix)
)

## Data list
data <- list(
  Lrel = tags$lenRel,
  Lrec = tags$lenRec,
  liberty = tags$libertyYears,
  Aoto=otoliths$age,
  Loto=otoliths$len,
  Lshort = 0,
  Llong = 60,
  Afix = 0.25
)

## Not estimated
map <- list(
  ## log_Linf = factor(NA),
  log_Lfix = factor(NA)
  ## log_sigma_1 = factor(NA),
  ## log_sigma_2 = factor(NA)
)

## Compile model
cp("boot/software/gcm.cpp", "model")
compile("model/gcm.cpp")
dyn.load(dynlib("model/gcm"))

## Run model
model <- MakeADFun(data, par, DLL="gcm", map=map, silent=TRUE)
fit <- nlminb(model$par, model$fn, model$gr,
              control=list(eval.max=1e4,iter.max=1e4))
report <- model$report()
sdreport <- sdreport(model, getReportCovariance=FALSE)

## Save results
save(model, fit, report, sdreport, file="model/gcm_1_results.RData")
