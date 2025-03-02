## Extract Richards results of interest, write TAF output tables

## Before: richards_results.RData (model)
## After:  richards_agepred.csv, richards_coefs.csv, richards_conv.csv,
##         richards_curve.csv, richards_residuals.csv,
##         richards_sigma.csv (output)

library(TAF)
library(RTMB)
library(fishgrowth)
source("utilities.R")  # ss

mkdir("output")

## Read model results
load("model/richards_results.RData")

## Extract coefficients
coefs <-
  as.data.frame(report[c("L1", "L2", "k", "b", "sigma_min", "sigma_max")])
L1 <- report$L1
L2 <- report$L2
k <- report$k
b <- report$b
sigma_min <- report$sigma_min
sigma_max <- report$sigma_max

## Calculate prediction sigma by length
L_min <- report$L_min
L_max <- report$L_max
sigma_slope <- (sigma_max - sigma_min) / (L_max - L_min)
sigma_intercept <- sigma_min - L_min * sigma_slope

## Calculate growth curve prediction band
t1 <- report$t1
t2 <- report$t2
Age <- seq(0, 4, 0.02)
Length <- richards_curve(Age, L1, L2, k, b, t1, t2)
curve <- data.frame(Age, Length)
curve$sigmaPred <- sigma_intercept + sigma_slope * curve$Length
curve$loPred <- curve$Length + qnorm(0.025) * curve$sigmaPred
curve$upPred <- curve$Length + qnorm(0.975) * curve$sigmaPred

## Model estimates
agepred <- data.frame(ageRel=report$age,
                      ageRec=report$age+report$liberty,
                      lenRel=report$Lrel,
                      lenRec=report$Lrec)
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
write.taf(agepred, "output/richards_agepred.csv")
write.taf(coefs, "output/richards_coefs.csv")
write.taf(conv, "output/richards_conv.csv")
write.taf(curve, "output/richards_curve.csv")
write.taf(otofit, "output/richards_otofit.csv")
write.taf(res, "output/richards_residuals.csv")
write.taf(sigma, "output/richards_sigma.csv")
