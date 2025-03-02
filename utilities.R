library(TAF)  # lim

################################################################################
## Basic plots of tag recapture data (lenRel, lenRec, libertyDays)

## Plot 1c.ii in Eveson et al. (2015)
plotLength <- function(dat, main, ...)
{
  if(missing(main))
    main <- paste0(deparse(substitute(dat)), " (", nrow(dat), ")")
  ylim <- lim(c(table(dat$lenRel), table(dat$lenRec)), 1.05)
  plot(table(dat$lenRec), col=NA, ylim=ylim,
       main=main, xlab="Length (cm)", ylab="Frequency", xaxt="n", ...)
  axis(1)
  dat$lenRel <- as.integer(round(dat$lenRel))
  dat$lenRec <- as.integer(round(dat$lenRec))
  points(table(dat$lenRel), lwd=3, col="blue", ...)
  points(table(dat$lenRec), lwd=3, col=rgb(0,1,0,alpha=0.7), ...)
  par("usr")
}

## Plot 1c.iii in Eveson et al. (2015)
plotLiberty <- function(dat, main, ...)
{
  if(missing(main))
    main <- paste("Avg", round(mean(dat$libertyDays)))
  hist(dat$libertyDays, col="white", xlab="Days at liberty", main=main, ...)
}

## Plot 2c in Eveson et al. (2015)
plotGrowth <- function(dat, from=-Inf, to=Inf, loess=TRUE, lm=FALSE, main, ...)
{
  dat <- dat[dat$lenAvg>=from & dat$lenAvg<=to,]
  x <- dat$lenAvg
  y <- dat$growth1
  if(missing(main))
    main <- paste("Avg", formatC(mean(y), format="f", digits=3))
  plot(y~x, main=main,
       xlab="Average length (cm)", ylab="Growth rate (cm/day)", ...)
  if(lm)
    abline(lm(y~x))
  if(loess)
    lines(loess(y~x), lwd=2, col="red")
}

################################################################################
## Generic loess plot

## From 'arni' package
plot.loess <- function(x, lwd=2, line.col="black", xlim=NULL, ylim=NULL,
                       xlab=NULL, ylab=NULL, ...)
{
  model <- x
  x   <- model$x
  y   <- model$y
  fit <- model$fit

  if(is.null(xlim))
    xlim <- range(x)
  if(is.null(ylim))
    ylim <- c(min(y,fit), max(y,fit))

  if(is.null(xlab))
    xlab <- as.character(terms(model))[3]
  if(is.null(ylab))
    ylab <- as.character(terms(model))[2]

  plot(x, y, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
  lines(x[order(x)], fit[order(x)], lwd=lwd, col=line.col)

  invisible(NULL)
}

## From 'arni' package
lines.loess <- function(x, ...)
{
  model <- x
  x <- model$x
  fit <- model$fit
  lines(x[order(x)], fit[order(x)], ...)
}

################################################################################
## Show model convergence statistics as a plot

plot.conv <- function(conv, x=c(-0.25,-0.15),
                      y=cumsum(c(0.9,rep(c(-0.05,-0.10),length=11))), cex=0.9)
{
  x <- rep(x, length=12)
  y <- rep(y, length=12)
  conv$iterations <- paste(c(conv$iterations, conv$evaluations.function,
                             conv$evaluations.gradient), collapse=", ")
  frame()
  text( x[1],  y[1], pos=4, xpd=TRUE, cex=cex, "objective")
  text( x[2],  y[2], pos=4, xpd=TRUE, cex=cex, round(conv$objective, 3))
  text( x[3],  y[3], pos=4, xpd=TRUE, cex=cex, "r2")
  text( x[4],  y[4], pos=4, xpd=TRUE, cex=cex, round(conv$r2, 3))
  text( x[5],  y[5], pos=4, xpd=TRUE, cex=cex, "convergence")
  text( x[6],  y[6], pos=4, xpd=TRUE, cex=cex, conv$convergence)
  text( x[7],  y[7], pos=4, xpd=TRUE, cex=cex, "iterations")
  text( x[8],  y[8], pos=4, xpd=TRUE, cex=cex, conv$iterations)
  text( x[9],  y[9], pos=4, xpd=TRUE, cex=cex, "message")
  text(x[10], y[10], pos=4, xpd=TRUE, cex=cex, conv$message)
  text(x[11], y[11], pos=4, xpd=TRUE, cex=cex, "pdHess")
  text(x[12], y[12], pos=4, xpd=TRUE, cex=cex, conv$pdHess)
}

ss <- function(x) sum((x-mean(x))^2)
