## Prepare GCM plots and tables for report

## Before: otoliths.csv (data), gcm_3_agepred.csv, gcm_3_coefs.csv,
##         gcm_3_conv.csv, gcm_3_curve.csv, gcm_3_residuals.csv (output)
## After:  gcm_3.pdf (report)

library(TAF)
library(areaplot)      # confplot
source("utilities.R")  # plot.conv

mkdir("report")

agepred <- read.taf("output/gcm_3_agepred.csv")
coefs <- read.taf("output/gcm_3_coefs.csv")
conv <- read.taf("output/gcm_3_conv.csv")
curve <- read.taf("output/gcm_3_curve.csv")
otoliths <- read.taf("data/otoliths.csv")
res <- read.taf("output/gcm_3_residuals.csv")

blue <- paste0(palette()[4], "b0")
red <- paste0(palette()[2], "b0")
black <- "#000000b0"

pdf("report/gcm_3.pdf", width=10, height=4)
par(mfrow=c(1,3))
plot(NA, xlim=c(0,4), ylim=c(20,85), xlab="Age (yr)", ylab="Length (cm)")
title(main=paste0("L0 = ", round(coefs$L0),
                  "    L4 = ", round(curve$Length[curve$Age == 4])))
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
