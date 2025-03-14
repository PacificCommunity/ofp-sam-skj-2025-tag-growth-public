## Compare growth curves

## Before: gcm_conv.csv, gcm_curve.csv, richards_conv.csv, richards_curve.csv,
##         schnute_3_conv.csv, schnute_3_curve.csv, vonbert_conv.csv,
##         vonbert_curve.csv (output)
## After:  comparison.csv, comparison.pdf (report)

library(TAF)

mkdir("report")

# Read convergence statistics
conv <- list()
conv$gcm <- read.taf("output/gcm_conv.csv")
conv$richards <- read.taf("output/richards_conv.csv")
conv$schnute_3 <- read.taf("output/schnute_3_conv.csv")
conv$vonbert <- read.taf("output/vonbert_conv.csv")

# Read growth curve
curve <- list()
curve$gcm <- read.taf("output/gcm_curve.csv")
curve$richards <- read.taf("output/richards_curve.csv")
curve$schnute_3 <- read.taf("output/schnute_3_curve.csv")
curve$vonbert <- read.taf("output/vonbert_curve.csv")

# Read coefficients
coefs <- list()
coefs$gcm <- read.taf("output/gcm_coefs.csv")
coefs$richards <- read.taf("output/richards_coefs.csv")
coefs$schnute_3 <- read.taf("output/schnute_3_coefs.csv")
coefs$vonbert <- read.taf("output/vonbert_coefs.csv")

# Plot growth curves
pdf("report/comparison.pdf", width=8, height=8)
plot(NA, xlim=c(0,4), ylim=c(20,85), xlab="Age (yr)", ylab="Length (cm)")
lines(Length~Age, curve$gcm, lwd=2, col=3)
lines(Length~Age, curve$richards, lwd=2, col=5)
lines(Length~Age, curve$schnute_3, lwd=2, col=6)
lines(Length~Age, curve$vonbert, lwd=2, col=7)
legend("bottomright", names(curve), lwd=3, col=1:7, bty="n", inset=0.02,
       y.intersp=1.25)
dev.off()

# Tabulate likelihood and parameter estimates
comparison <- data.frame(model=names(conv), coefs=c(4,4,3,3),
                         nll=round(sapply(conv, `[[`, "objective"), 3),
                         hessian=sapply(conv, `[[`, "pdHess"), row.names=NULL)

comparison$est[comparison$model=="gcm"] <-
  paste0("L0=", round(coefs$gcm$L0, 3),
         ", rmax=", round(coefs$gcm$rmax, 3),
         ", k=", round(coefs$gcm$k, 3),
         ", t50=", round(coefs$gcm$t50, 3))
comparison$sigma[comparison$model=="gcm"] <-
  paste(formatC(coefs$gcm$sigma_min, format="f", digits=3),
        formatC(coefs$gcm$sigma_max, format="f", digits=3), sep=", ")

comparison$est[comparison$model=="richards"] <-
  paste0("L1=", round(coefs$richards$L1, 3),
         ", L2=", round(coefs$richards$L2, 3),
         ", k=", formatC(coefs$richards$k, format="f", digits=3),
         ", b=", round(coefs$richards$b, 3))
comparison$sigma[comparison$model=="richards"] <-
  paste(formatC(coefs$richards$sigma_min, format="f", digits=3),
        formatC(coefs$richards$sigma_max, format="f", digits=3), sep=", ")

comparison$est[comparison$model=="schnute_3"] <-
  paste0("L1=", round(coefs$schnute_3$L1, 3),
         ", L2=", round(coefs$schnute_3$L2, 3),
         ", b=", round(coefs$schnute_3$b, 3))
comparison$sigma[comparison$model=="schnute_3"] <-
  paste(formatC(coefs$schnute_3$sigma_min, format="f", digits=3),
        formatC(coefs$schnute_3$sigma_max, format="f", digits=3), sep=", ")

comparison$est[comparison$model=="vonbert"] <-
  paste0("L1=", round(coefs$vonbert$L1, 3),
         ", L2=", formatC(coefs$vonbert$L2, format="f", digits=3),
         ", k=", round(coefs$vonbert$k, 3))
comparison$sigma[comparison$model=="vonbert"] <-
  paste(formatC(coefs$vonbert$sigma_min, format="f", digits=3),
        formatC(coefs$vonbert$sigma_max, format="f", digits=3), sep=", ")

write.taf(comparison, dir="report", quote=TRUE)
