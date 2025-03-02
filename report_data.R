## Prepare data plots and tables for report

## Before: otoliths.csv, tags_jp.csv, tags_shi.csv (data)
## After:  data_jp_shi.pdf, data_tags_inc.pdf, data_tags_oto.pdf (report)

library(TAF)
source("utilities.R")  # plotGrowth, plotLength, plotLiberty

mkdir("report")

otoliths <- read.taf("data/otoliths.csv")
jp <- read.taf("data/tags_jp.csv")
shi <- read.taf("data/tags_shi.csv")
tags <- rbind(jp[intersect(names(jp), names(shi))],
              shi[intersect(names(jp), names(shi))])

## JPTP and PTTP-shimizu
pdf("report/data_jp_shi.pdf", width=12, height=6)
par(mfrow=c(2,4))
plotLength(jp, xlim=c(30,80))
plotLiberty(jp, nclass=20, xlim=c(0,1100))
plotGrowth(jp, to=70, xlim=c(30,70), ylim=c(0,0.2))
plotGrowth(jp, from=40, to=60, xlim=c(40,60), ylim=c(0,0.2), lm=TRUE)
plotLength(shi, xlim=c(30,80))
plotLiberty(shi, nclass=6, xlim=c(0,1100))
plotGrowth(shi, to=70, xlim=c(30,70), ylim=c(0,0.2))
plotGrowth(shi, from=40, to=60, xlim=c(40,60), ylim=c(0,0.2), lm=TRUE)
dev.off()

## Tags (combined JPTP and PTTP-shimizu) and otoliths
pdf("report/data_tags_oto.pdf", width=12, height=6)
par(mfrow=c(2,4))
## combined tags
plotLength(tags, xlim=c(30,80))
plotLiberty(tags, nclass=20, xlim=c(0,1100))
plotGrowth(tags, to=70, xlim=c(30,70), ylim=c(0,0.2))
plotGrowth(tags, from=40, to=60, xlim=c(40,60), ylim=c(0,0.2), lm=TRUE)
## otoliths
hist(otoliths$len, xlab="Length (cm)", xlim=c(30,80), col="white",
     main=paste0("otoliths (", nrow(otoliths), ")"))
hist(otoliths$age, xlab="Age (yr)", col="white",
     main=paste("Avg", round(mean(otoliths$age),1)))
plot(len~age, otoliths, xlab="Age (yr)", ylab="Length (cm)")
abline(lm(len~age, otoliths))
avg.oto <- coef(lm(len~age, otoliths))[[2]] / 365
main <- paste("Avg", formatC(mean(avg.oto), format="f", digits=3))
title(main=main)
plot(len~age, otoliths, subset=len>=40 & len<=60, ylim=c(40,60),
     xlab="Age (yr)", ylab="Length (cm)")
abline(lm(len~age, otoliths, subset=len>=40 & len<=60))
avg.oto <- coef(lm(len~age, otoliths, subset=len>=40 & len<=60))[[2]] / 365
main <- paste("Avg", formatC(mean(avg.oto), format="f", digits=3))
title(main=main)
dev.off()

pdf("report/data_tags_inc.pdf", width=10, height=6)
xlim <- lim(c(jp$libertyDays, shi$libertyDays), 1.05)
ylim <- lim(c(jp$lenDiff, shi$lenDiff), 1.05)
plot(NA, xlim=xlim, ylim=ylim,
     xlab="Days at liberty", ylab="Length increment (cm)")
points(lenDiff~libertyDays, jp, col="sienna")
points(lenDiff~libertyDays, shi, col=3)
legend("bottomright", c("JPTP","PTTP (Shimizu)"), pch=1, col=c("sienna",3),
       bty="n", inset=0.02, y.intersp=1.1)
dev.off()
