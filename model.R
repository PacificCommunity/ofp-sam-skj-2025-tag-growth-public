## Run analysis, write model results

library(TAF)

sourceTAF("model_rate.R")
sourceTAF("model_gcm_1.R"); detach("package:TMB")
sourceTAF("model_gcm_2.R")
sourceTAF("model_gcm_3.R"); detach("package:fishgrowth"); detach("package:RTMB")
sourceTAF("model_gcm_4.R"); detach("package:TMB")
sourceTAF("model_richards.R")
sourceTAF("model_schnute_3.R")
sourceTAF("model_vonbert.R")
