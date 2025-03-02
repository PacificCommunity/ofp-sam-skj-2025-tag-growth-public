## Prepare data plots and tables for report

## Before: otoliths.csv (data), rate.csv (output)
## After:  rate.pdf (report)

library(TAF)

mkdir("report")

otoliths <- read.taf("data/otoliths.csv")
rate <- read.taf("output/rate.csv")

# Plot otoliths
pdf("report/rate_slope.pdf")
plot(otoliths, xlim=c(0,1.5), ylim=c(0,80),
     main="Growth rate from otoliths and tags", pch=16)
fm <- lm(len~age, otoliths)
abline(fm)
oto.rate <- 10 * coef(fm)[[2]] / 365  # mm/day
# growth rate slope from tags
tags <- rate[rate$Data == "Tags",]
abline(a=37, b=365*mean(tags$Rate)/10, lty=2)
legend("bottomright", c("Otoliths","Tags"), lwd=1.5, lty=1:2, bty="n",
       inset=0.02, y.intersp=1.25)
dev.off()

# Boxplot
pdf("report/rate_boxplot.pdf", width=6)
boxplot(tags$Rate, main="Growth rate from otoliths and tags",
        ylab="daily growth rate (mm/day)")
abline(h=oto.rate, lwd=2, col=2)
opar <- par(lend="butt")
legend("bottomright", c("Otoliths","Tags"), col=c(2,"gray"), lwd=c(2.5, 12),
       bty="n", inset=0.02, y.intersp=1.5)
par(opar)
box()
dev.off()
