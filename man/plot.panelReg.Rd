\name{plot.panelReg}
\alias{plot.panelReg}
\title{
  Plot an panelReg Object.
}
\description{
  Plot the estimated coefficient function(s) from a fit
  of Aalen's additive regression model.
}
\usage{
\method{plot}{panelReg}(x, se=TRUE, ...)
}
\arguments{
  \item{x}{
    The result of a call to the \code{panelReg} function
  }
  \item{se}{
    If \code{TRUE}, standard error bands are included on the plot.
  }
  \item{\dots }{
    Other graphical parameters such as line type, color, or axis labels.
  }
}
\examples{
# Plot the fit of bladder tumor data set
formula <- PanelSurv(id, time, count) ~ num + size + treatment

fit <- panelReg(formula, data=blaTum, method="AEE", se="Sandwich")
fit

plot(fit, se=TRUE)
}
\seealso{
  \code{\link{panelReg}}, \code{\link{panelReg.object}}.
}