\name{bladTumor}
\docType{data}
\alias{bladTumor}
\title{Bladder Tumors Cancer Recurrences}
\usage{
bladTumor
}
\description{A data frame contains data on recurrences of bladder
  cancer, used by many people to demonstrate methodology for recurrent
  event modelling. The data was obtained by courtesy of Ying Zhang,
  containing records of 118 patients from three treatment arms: 48 are
  from the placebo arm, 37 are from the thiotepa arm, and the rest 33
  are from the pyridoxine arm.
}
\format{
  \describe{
    \item{\code{subject}}{: patient id}
    \item{\code{time}}{: observation time}
    \item{\code{count}}{: cumulative number of tumors}
    \item{\code{count2}}{: number of new tumors since last observation time}
    \item{\code{number}}{: initial number of tumours (8=8 or more)}
    \item{\code{size}}{: size (cm) of largest initial tumour}
    \item{\code{pyridoxine}}{: dummy variable for pyridoxine arm}
    \item{\code{thiotepa}}{: dummy variable for thiotepa arm}
  }
}
\references{
  Byar, D. P. (1980). The Veterans administration study of
  chemoprophylaxis for recurrent stage I bladder tumors: comparisons
  of placebo, pyridoxine, and topical thiotepa. \emph{Bladder Tumors and
  Other Topics in Urological Oncology}, pp. 363--370. New York: Plenum.

  Wellner, J. A. and Zhang, Y. (2007) Two likelihood-based
  semiparametric estimation methods for panel count data with
  covariates. \emph{Annals of Statistics}, \bold{35} 2106--2142.

  Lu, M., Zhang, Y. and Huang, J. (2009) Semiparametric estimation methods for
  panel count data using monotone B-Splines. \emph{Journal of the American
  Statistical Association} \bold{104} 1060--1070.
}
\note{
  To our surprise, the two-treatment (placebo and thiotepa) subset of
  the full version \code{bladTumor} do not match the two-treatment
  version \code{blaTum}.
}
\seealso{
  \code{\link{blaTum}}
}
\examples{
# Load data set
data(bladTumor)

# Plot bladder tumor data
p <- plot(with(bladTumor, PanelSurv(subject, time, count2)))
print(p)

\dontrun{
# Fit a semiparametric regression model
formula <- PanelSurv(subject, time, count2) ~ number + size + pyridoxine + thiotepa

panelReg(formula, data=bladTumor, method="MPL", se="Bootstrap", control=list(R=20))

panelReg(formula, data=bladTumor, method="MLs", se="Bootstrap", control=list(R=20))
}
}
\keyword{survival}
\keyword{panel count}
