## Prepare von Bertalanffy plots and tables for report

## Before: otoliths.csv (data), vonbert_agepred.csv,
##         vonbert_coefs.csv, vonbert_conv.csv,
##         vonbert_curve.csv, vonbert_residuals.csv (output)
## After:  vonbert.pdf (report)

library(TAF)
library(areaplot)      # confplot
source("utilities.R")  # plot.conv

mkdir("report")

agepred <- read.taf("output/vonbert_agepred.csv")
coefs <- read.taf("output/vonbert_coefs.csv")
conv <- read.taf("output/vonbert_conv.csv")
curve <- read.taf("output/vonbert_curve.csv")
otoliths <- read.taf("data/otoliths.csv")
res <- read.taf("output/vonbert_residuals.csv")

blue <- paste0(palette()[4], "b0")
red <- paste0(palette()[2], "b0")
black <- "#000000b0"

pdf("report/vonbert.pdf", width=10, height=4)
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
title(main=paste("k =", round(coefs$k, 2)))
points(Residual~Age, res, subset=Component=="TagRel", col=blue)
points(Residual~Age, res, subset=Component=="TagRec", col=red)
points(Residual~Age, res, subset=Component=="Otoliths", col=black)
abline(h=0)

plot.conv(conv)
dev.off()
