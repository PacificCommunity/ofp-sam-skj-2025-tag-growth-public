## Preprocess otolith data, write TAF data tables

## Before: Leroy_oto_tag_estimates.csv (boot/data),
## After:  otoliths.csv (data)

library(TAF)

mkdir("data")

## Read original data
otoliths <- read.taf("boot/data/Leroy_oto_tag_estimates.csv")

## Filter rows and columns
otoliths <- otoliths[otoliths$source == "otolith", c("age_yrs","FL_cm")]

## Rename columns
names(otoliths) <- c("age", "len")

## Filter
otoliths <- otoliths[otoliths$len <= 45,]

## Write table
write.taf(otoliths, dir="data")
