## Extract results of interest, write TAF output tables

library(TAF)

mkdir("output")

cp("model/rate.csv", "output")

sourceTAF("output_gcm_1.R")
sourceTAF("output_gcm_2.R")
sourceTAF("output_gcm_3.R")
sourceTAF("output_gcm_4.R")
sourceTAF("output_richards.R")
sourceTAF("output_schnute_3.R")
sourceTAF("output_vonbert.R")
