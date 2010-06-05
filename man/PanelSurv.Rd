\name{PanelSurv}
\alias{PanelSurv}
\alias{is.PanelSurv}
\alias{plot.PanelSurv}
\title{Create a PanelSurv Object}
\description{
  Create a panel count survival object, usually used as a response
  variable in a model formula.
}
\usage{
PanelSurv(ID, time, count)
\method{is}{PanelSurv}(x)
\method{plot}{PanelSurv}(x, ...)
}
\arguments{
  \item{ID}{Observation subject's ID.}
  \item{time}{Observation time.}
  \item{count}{Number of new events since last observation time.}
  \item{x}{An \code{PanelSurv} object.}
  \item{\dots }{
    Other graphical parameters such as line type, color, or axis labels.
  }
}
%\details{
%  ~~ to add ~~
%}
\value{
  An object of S3 class \code{"PanelSurv"}.
  \item{psDF}{a data frame, part of original input data frame with variable "ID", "time"
    and "count".}
  \item{timeGrid}{ordered distinct observation times in the set of all
    observation times.}
  \item{panelMatrix}{a matrix representation of panel count data, one
    row per subject, one column per time point in \code{"timeGrid"}.}

  In the case of \code{is.PanelSurv}, a logical value \code{TRUE} if
  \code{x} inherits from class \code{"PanelSurv"}, otherwise an \code{FALSE}.

  In the case of \code{plot.PanelSurv}, a tile plot of
  \code{panelMatrix} produced by package \code{ggplot2} with color
  indicating number of counts since last observation time.
}
% \references{ ~put references to the literature/web site here ~ }

\seealso{\code{\link{panelReg}}}
\examples{
data(blaTum)

response <- with(blaTum, PanelSurv(id, time, count))
is.PanelSurv(response)
plot(response)
}
\keyword{survival}
\keyword{panel count}