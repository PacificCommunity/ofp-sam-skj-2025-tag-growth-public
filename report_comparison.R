## Compare growth curves

## Before: gcm_1_conv.csv, gcm_1_curve.csv, gcm_2_conv.csv, gcm_2_curve.csv,
##         gcm_3_conv.csv, gcm_3_curve.csv, gcm_4_conv.csv, gcm_4_curve.csv,
##         richards_conv.csv, richards_curve.csv, schnute_3_conv.csv,
##         schnute_3_curve.csv, vonbert_conv.csv, vonbert_curve.csv (output)
## After:  comparison.csv, comparison.pdf (report)

library(TAF)

mkdir("report")

# Read convergence statistics
conv <- list()
conv$gcm_1 <- read.taf("output/gcm_1_conv.csv")
conv$gcm_2 <- read.taf("output/gcm_2_conv.csv")
conv$gcm_3 <- read.taf("output/gcm_3_conv.csv")
conv$gcm_4 <- read.taf("output/gcm_4_conv.csv")
conv$richards <- read.taf("output/richards_conv.csv")
conv$schnute_3 <- read.taf("output/schnute_3_conv.csv")
conv$vonbert <- read.taf("output/vonbert_conv.csv")

# Read growth curve
curve <- list()
curve$gcm_1 <- read.taf("output/gcm_1_curve.csv")
curve$gcm_2 <- read.taf("output/gcm_2_curve.csv")
curve$gcm_3 <- read.taf("output/gcm_3_curve.csv")
curve$gcm_4 <- read.taf("output/gcm_4_curve.csv")
curve$richards <- read.taf("output/richards_curve.csv")
curve$schnute_3 <- read.taf("output/schnute_3_curve.csv")
curve$vonbert <- read.taf("output/vonbert_curve.csv")

# Read coefficients
coefs <- list()
coefs$gcm_1 <- read.taf("output/gcm_1_coefs.csv")
coefs$gcm_2 <- read.taf("output/gcm_2_coefs.csv")
coefs$gcm_3 <- read.taf("output/gcm_3_coefs.csv")
coefs$gcm_4 <- read.taf("output/gcm_4_coefs.csv")
coefs$richards <- read.taf("output/richards_coefs.csv")
coefs$schnute_3 <- read.taf("output/schnute_3_coefs.csv")
coefs$vonbert <- read.taf("output/vonbert_coefs.csv")

# Plot growth curves
pdf("report/comparison.pdf", width=8, height=8)
plot(NA, xlim=c(0,4), ylim=c(20,85), xlab="Age (yr)", ylab="Length (cm)")
lines(Length~Age, curve$gcm_1, lwd=2, col=1)
lines(Length~Age, curve$gcm_2, lwd=2, col=2)
lines(Length~Age, curve$gcm_3, lwd=2, col=3)
lines(Length~Age, curve$gcm_4, lwd=2, col=4)
lines(Length~Age, curve$richards, lwd=2, col=5)
lines(Length~Age, curve$schnute_3, lwd=2, col=6)
lines(Length~Age, curve$vonbert, lwd=2, col=7)
legend("bottomright", names(curve), lwd=3, col=1:7, bty="n", inset=0.02,
       y.intersp=1.25)
dev.off()

# Tabulate likelihood and parameter estimates
comparison <- data.frame(model=names(conv), coefs=c(4,4,4,4,4,3,3),
                         nll=round(sapply(conv, `[[`, "objective"), 3),
                         hessian=sapply(conv, `[[`, "pdHess"), row.names=NULL)

comparison$est[comparison$model=="gcm_1"] <-
  paste0("Linf=", round(coefs$gcm_1$Linf, 3),
         ", k=", round(coefs$gcm_1$k, 3),
         ", rmax=", round(coefs$gcm_1$rmax, 3),
         ", Lfix=34*")
comparison$sigma[comparison$model=="gcm_1"] <-
  paste(formatC(coefs$gcm_1$sigma_1, format="f", digits=3),
        formatC(coefs$gcm_1$sigma_2, format="f", digits=3), sep=", ")

comparison$est[comparison$model=="gcm_2"] <-
  paste0("L0=", round(coefs$gcm_2$L0, 3),
         ", rmax=", round(coefs$gcm_2$rmax, 3),
         ", k=", round(coefs$gcm_2$k, 3),
         ", t50=", round(coefs$gcm_2$t50, 3))
comparison$sigma[comparison$model=="gcm_2"] <-
  paste(formatC(coefs$gcm_2$sigma_min, format="f", digits=3),
        formatC(coefs$gcm_2$sigma_max, format="f", digits=3), sep=", ")

comparison$est[comparison$model=="gcm_3"] <-
  paste0("L0=", round(coefs$gcm_3$L0, 3),
         ", rmax=", round(coefs$gcm_3$rmax, 3),
         ", k=", round(coefs$gcm_3$k, 3),
         ", t50=", round(coefs$gcm_3$t50, 3))
comparison$sigma[comparison$model=="gcm_3"] <-
  paste(formatC(coefs$gcm_3$sigma_min, format="f", digits=3),
        formatC(coefs$gcm_3$sigma_max, format="f", digits=3), sep=", ")

comparison$est[comparison$model=="gcm_4"] <-
  paste0("Linf=", round(coefs$gcm_4$Linf, 3),
         ", k=", round(coefs$gcm_4$k, 3),
         ", rmax=", round(coefs$gcm_4$rmax, 3),
         ", Lfix=34*")
comparison$sigma[comparison$model=="gcm_4"] <-
  paste(formatC(coefs$gcm_4$sigma_1, format="f", digits=3),
        formatC(coefs$gcm_4$sigma_2, format="f", digits=3), sep=", ")

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
