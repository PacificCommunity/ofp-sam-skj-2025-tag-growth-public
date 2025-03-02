## Extract GCM results of interest, write TAF output tables

## Before: gcm_4_results.RData (model)
## After:  gcm_4_agepred.csv, gcm_4_coefs.csv, gcm_4_conv.csv, gcm_4_curve.csv,
##         gcm_4_otofit.csv, gcm_4_residuals.csv, gcm_4_sigma.csv (output)

library(TAF)
source("utilities.R")  # ss

mkdir("output")

## Read model results
load("model/gcm_4_results.RData")

## Extract coefficients
coefs <- as.data.frame(report[c("Linf", "k", "rmax", "sigma_1", "sigma_2")])
Linf <- report$Linf
k <- report$k
rmax <- report$rmax
sigma.1 <- report$sigma_1
sigma.2 <- report$sigma_2
L0 <- report$Lfix - rmax * report$Afix
t50 <- log(exp(k * (Linf-L0) / rmax) - 1) / k

## Calculate prediction sigma by length
Lshort <- report$Lshort
Llong <- report$Llong
sigma.slope <- (sigma.2 - sigma.1) / (Llong - Lshort)
sigma.intercept <- sigma.1 - Lshort * sigma.slope

## Calculate growth curve prediction bands
t1 <- report$t1
t2 <- report$t2
Age <- seq(0, 4, 0.02)
Length <- L0 + rmax * ((log(exp(-k*t50)+1) - log(exp(k*(Age-t50))+1)) / k + Age)
curve <- data.frame(Age, Length)
curve$sigmaPred <- sigma.intercept + sigma.slope * curve$Length
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
sigmaPred <- sigma.intercept + sigma.slope * res$Length_hat
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
midpoint <- mean(c(report$Lshort, report$Llong))
empiricalYoung <- sd(res$Residual[res$Length_hat < midpoint])
empiricalOld <- sd(res$Residual[res$Length_hat >= midpoint])
sigma <- data.frame(empiricalYoung, empiricalOld, coverage)

## Export tables
write.taf(agepred, "output/gcm_4_agepred.csv")
write.taf(coefs, "output/gcm_4_coefs.csv")
write.taf(conv, "output/gcm_4_conv.csv")
write.taf(curve, "output/gcm_4_curve.csv")
write.taf(otofit, "output/gcm_4_otofit.csv")
write.taf(res, "output/gcm_4_residuals.csv")
write.taf(sigma, "output/gcm_4_sigma.csv")
