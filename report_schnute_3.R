## Prepare Schnute case 3 plots and tables for report

## Before: otoliths.csv (data), schnute_3_agepred.csv, schnute_3_coefs.csv,
##         schnute_3_conv.csv, schnute_3_curve.csv,
##         schnute_3_residuals.csv (output)
## After:  schnute_3.pdf (report)

library(TAF)
library(areaplot)      # confplot
source("utilities.R")  # plot.conv

mkdir("report")

agepred <- read.taf("output/schnute_3_agepred.csv")
coefs <- read.taf("output/schnute_3_coefs.csv")
conv <- read.taf("output/schnute_3_conv.csv")
curve <- read.taf("output/schnute_3_curve.csv")
otoliths <- read.taf("data/otoliths.csv")
res <- read.taf("output/schnute_3_residuals.csv")

blue <- paste0(palette()[4], "b0")
red <- paste0(palette()[2], "b0")
black <- "#000000b0"

pdf("report/schnute_3.pdf", width=10, height=4)
par(mfrow=c(1,3))
plot(NA, xlim=c(0,4), ylim=c(20,85), xlab="Age (yr)", ylab="Length (cm)")
title(main=paste0("L0 = ", round(coefs$L1), "    L4 = ", round(coefs$L2)))
confplot(curve[c("Age", "loPred", "upPred")], add=TRUE)
points(lenRel~ageRel, data=agepred, col=blue)
points(lenRec~ageRec, data=agepred, col=red)
lines(Length~Age, curve)
points(otoliths, col=black)

plot(NA, xlim=c(0,4), ylim=c(-16,16),
     xlab="Age (yr)", ylab="Length residual (cm)")
title(main=paste("b =", formatC(round(coefs$b, 2), format="f", digits=2)))
points(Residual~Age, res, subset=Component=="TagRel", col=blue)
points(Residual~Age, res, subset=Component=="TagRec", col=red)
points(Residual~Age, res, subset=Component=="Otoliths", col=black)
abline(h=0)

plot.conv(conv)
dev.off()
