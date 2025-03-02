## Analyze daily growth rate data, write model results

## Before: otoliths.csv, tags_jp.csv, tags_shi.csv (data)
## After:  rate.csv (model)

library(TAF)

mkdir("model")

# Read data
otoliths <- read.taf("data/otoliths.csv")
tags.jp <- read.taf("data/tags_jp.csv")
tags.shi <- read.taf("data/tags_shi.csv")

# Calculate daily growth rate (mm/day)
fm <- lm(len~age, otoliths)
otoliths <- 10 * coef(fm)[[2]] / 365
tags <- c(tags.jp$growth1[tags.jp$lenAvg<=45],
          tags.shi$growth1[tags.shi$lenAvg<=45])
tags <- 10 * tags

# Construct data frame
otoliths.df <- data.frame(Data="Otoliths", Subset="All", Rate=otoliths)
tags.df <- data.frame(Data="Tags",
                      Subset=rep(c("JP", "Shi"),
                                 c(sum(tags.jp$lenAvg<=45),
                                   sum(tags.shi$lenAvg<=45))), Rate=tags)
rate <- rbind(otoliths.df, tags.df)

# Write TAF table
write.taf(rate, dir="model")
