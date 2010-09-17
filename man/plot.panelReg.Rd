\name{plot.panelReg}
\alias{plot.panelReg}
\title{
  Plot a panelReg Object.
}
\description{
  Plot the estimated baseline mean function. If \code{"se"} option of
  \code{panelReg} is not \code{"NULL"}, 95 percent point wise confidence
  interval is also ploted.
}
\usage{
\method{plot}{panelReg}(x, ...)
}
\arguments{
  \item{x}{
    The result of a call to the \code{panelReg} function.
  }
  \item{\dots}{
    Other graphical parameters such as line type, color, or axis labels.
  }
}
\examples{
data(blaTum)

# Plot the fit of bladder tumor data set
formula <- PanelSurv(id, time, count) ~ num + size + treatment

fit1 <- panelReg(formula, data=blaTum, method="AEE", se="Sandwich")
plot(fit1)

fit2 <- panelReg(formula, data=blaTum, method="MLs", se="NULL")
plot(fit2)
}
\seealso{
  \code{\link{panelReg}}, \code{\link{panelReg.object}}.
}