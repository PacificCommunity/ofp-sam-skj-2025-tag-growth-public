## Extract GCM results of interest, write TAF output tables

## Before: gcm_3_results.RData (model)
## After:  gcm_3_agepred.csv, gcm_3_coefs.csv, gcm_3_conv.csv, gcm_3_curve.csv,
##         gcm_3_residuals.csv, gcm_3_sigma.csv (output)

library(TAF)
library(RTMB)
library(fishgrowth)
source("utilities.R")  # ss

mkdir("output")

## Read model results
load("model/gcm_3_results.RData")

## Extract coefficients
coefs <-
  as.data.frame(report[c("L0", "rmax", "k", "t50", "sigma_min", "sigma_max")])
L0 <- report$L0
rmax <- report$rmax
k <- report$k
t50 <- report$t50
sigma_min <- report$sigma_min
sigma_max <- report$sigma_max

## Calculate prediction sigma by length
L_min <- report$L_min
L_max <- report$L_max
sigma_slope <- (sigma_max - sigma_min) / (L_max - L_min)
sigma_intercept <- sigma_min - L_min * sigma_slope

## Calculate growth curve prediction band
Age <- seq(0, 4, 0.02)
Length <- gcm_curve(Age, L0, rmax, k, t50)
curve <- data.frame(Age, Length)
curve$sigmaPred <- sigma_intercept + sigma_slope * curve$Length
curve$loPred <- curve$Length + qnorm(0.025) * curve$sigmaPred
curve$upPred <- curve$Length + qnorm(0.975) * curve$sigmaPred

## Model estimates
agepred <- data.frame(ageRel=report$age,
                      ageRec=report$age+report$liberty,
                      lenRel=report$Lrel, lenRec=report$Lrec)
otofit <- data.frame(Age=report$Aoto,
                     Length=report$Loto,
                     Length_hat=report$Loto_hat,
                     Residual=report$Loto-report$Loto_hat)

## Residuals
res <- data.frame(Component=c(rep("TagRel", nrow(agepred)),
                              rep("TagRec", nrow(agepred)),
                              rep("Otoliths", nrow(otofit))),
                  Age=c(agepred$ageRel, agepred$ageRec, otofit$Age),
                  Length=c(agepred$lenRel, agepred$lenRec, otofit$Length),
                  Length_hat=c(report$Lrel_hat, report$Lrec_hat,
                               otofit$Length_hat))
res$Residual <- res$Length - res$Length_hat

## Prediction band coverage
sigmaPred <- sigma_intercept + sigma_slope * res$Length_hat
loPred <- res$Length_hat + qnorm(0.025) * sigmaPred
upPred <- res$Length_hat + qnorm(0.975) * sigmaPred
inside <- res$Length >= loPred & res$Length <= upPred
res$Inside95band <- inside
coverage <- 100 * round(sum(inside)/length(inside), 3)

## Convergence and goodness of fit
conv <- as.data.frame(as.list(unlist(fit[-1])))
conv$pdHess <- sdreport$pdHess
rel <- res$Length[res$Component=="TagRel"]
rec <- res$Length[res$Component=="TagRec"]
oto <- res$Length[res$Component=="Otoliths"]
total <- ss(rel) + ss(rec) + ss(oto)
remaining <- sum(res$Residual^2)
explained <- total - remaining
conv$r2 <- explained / total

## Empirical sigma
midpoint <- mean(c(report$L_min, report$L_max))
empiricalYoung <- sd(res$Residual[res$Length_hat < midpoint])
empiricalOld <- sd(res$Residual[res$Length_hat >= midpoint])
sigma <- data.frame(empiricalYoung, empiricalOld, coverage)

## Export tables
write.taf(agepred, "output/gcm_3_agepred.csv")
write.taf(coefs, "output/gcm_3_coefs.csv")
write.taf(conv, "output/gcm_3_conv.csv")
write.taf(curve, "output/gcm_3_curve.csv")
write.taf(otofit, "output/gcm_3_otofit.csv")
write.taf(res, "output/gcm_3_residuals.csv")
write.taf(sigma, "output/gcm_3_sigma.csv")
