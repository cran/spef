# Internal functions :
#   doPanelFitPlot.stepfun
#   doPanelFitPlot.stepfun.stepfun
#   doPanelFitPlot.isplineFun
# Export functions :
#   plot.panelReg
#   print.panelReg

##############################################################################
doPanelFitPlot.stepfun <- function(baseline, timeGrid, baselineSE, ...) {
    plot(baseline, do.points=FALSE, xlab="Time",
         ylab=expression(hat(Lambda[0])(t)),
         main="Cumulative Baseline Mean", ...)
}

##############################################################################
doPanelFitPlot.stepfun.se <- function(baseline, timeGrid, baselineSE, ...) {
    y <- baseline(timeGrid)
    lowFun <- stepfun(timeGrid, c(0, y * exp(- 1.96 * baselineSE / y)))
    highFun <- stepfun(timeGrid, c(0, y * exp(1.96 * baselineSE / y)))

    plot(highFun, do.points=FALSE, lty=2, xlab="Time",
         ylab=expression(hat(Lambda[0])(t)), main="Cumulative Baseline Mean")
    plot(baseline, do.points=FALSE, add=TRUE, ...)
    plot(lowFun, do.points=FALSE, lty=2, add=TRUE)
}

##############################################################################
doPanelFitPlot.isplineFun <- function(baseline, timeGrid, baselineSE, ...) {
    plot(baseline, xlab="Time", ylab=expression(hat(Lambda[0])(t)),
         main="Cumulative Baseline Mean (I-Spline)", ...)
}

doPanelFitPlot.isplineFun.se <- function(baseline, timeGrid, baselineSE, ...) {
    y <- baseline(timeGrid)
    low <- y * exp(- 1.96 * baselineSE / y)
    high <- y * exp(1.96 * baselineSE / y)

    plot(baseline, xlab="Time", ylab=expression(hat(Lambda[0])(t)),
         main="Cumulative Baseline Mean (I-Spline)",
         ylim=c(0, 1.05 * max(high)), ...)
    points(timeGrid, high, type="l", lty=2)
    points(timeGrid, low, type="l", lty=2)
}

##############################################################################
# Method dispatch
##############################################################################
setGeneric("doPanelFitPlot",
           function(baseline, baselineSE, ...) {
               standardGeneric("doPanelFitPlot")
           })

setOldClass(c("stepfun", "function"))

setMethod("doPanelFitPlot",
          signature(baseline="stepfun", baselineSE="NULL"),
          doPanelFitPlot.stepfun)

setMethod("doPanelFitPlot",
          signature(baseline="stepfun", baselineSE="numeric"),
          doPanelFitPlot.stepfun.se)

setOldClass(c("isplineFun", "function"))

setMethod("doPanelFitPlot",
          signature(baseline="isplineFun", baselineSE="NULL"),
          doPanelFitPlot.isplineFun)

setMethod("doPanelFitPlot",
          signature(baseline="isplineFun", baselineSE="numeric"),
          doPanelFitPlot.isplineFun.se)


##############################################################################
# Plot a PanelReg object
##############################################################################
plot.panelReg <- function(x, ...) {
    doPanelFitPlot(baseline=x$baseline,
                   timeGrid=x$timeGrid,
                   baselineSE=x$baselineSE, ...)
}

##############################################################################
# Print a panelReg object
##############################################################################
print.panelReg <- function(x, digits=max(options()$digits - 4, 3), ...) {
    savedig <- options(digits = digits)
    on.exit(options(savedig))

    coef <- x$beta
    if (is.null(x$betaSE))
        se <- rep(NA, length(coef))
    else
        se <- x$betaSE

    tmp <- data.frame(coef, exp(coef), se, coef/se,
                      signif(1 - pchisq((coef/ se)^2, 1), digits -1))

    dimnames(tmp) <- list(names(coef),
                          c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)"))

    # Print results
    cat("\n")
    cat("Call:\n")
    dput(x$call)
    cat("\n")
    printCoefmat(tmp)
    cat("\n")
    invisible()
}
