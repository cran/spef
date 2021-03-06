% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spef_data.R
\docType{data}
\name{bladTumor}
\alias{bladTumor}
\title{Bladder Tumors Cancer Recurrences}
\format{A \code{data.frame} contains the following columns:
\describe{
  \item{\code{subject}}{patient id}
  \item{\code{time}}{observation time}
  \item{\code{count}}{cumulative number of tumors}
  \item{\code{count2}}{number of new tumors since last observation time}
  \item{\code{number}}{initial number of tumors (8=8 or more)}
  \item{\code{size}}{size (cm) of largest initial tumor}
  \item{\code{pyridoxine}}{dummy variable for pyridoxine arm}
  \item{\code{thiotepa}}{dummy variable for thiotepa arm}
}}
\usage{
data(bladTumor)
}
\description{
A data frame contains data on recurrences of bladder cancer,
used by many people to demonstrate methodology for recurrent event modeling.
The data was obtained by courtesy of Ying Zhang, containing records of 118 patients
from three treatment arms: 48 are from the placebo arm, 37 are from the thiotepa arm,
and the rest 33 are from the pyridoxine arm.
}
\note{
To our surprise, the two-treatment (placebo and thiotepa) subset of
  the full version \code{bladTumor} do not match the two-treatment
  version \code{blaTum}.
}
\examples{
data(bladTumor)
## Plot bladder tumor data
p <- plot(with(bladTumor, PanelSurv(subject, time, count2)))
print(p)
}
\references{
Byar, D. P. (1980). The Veterans administration study of chemoprophylaxis
for recurrent stage I bladder tumors: Comparisons of placebo, pyridoxine, and topical thiotepa.
\emph{Bladder Tumors and Other Topics in Urological Oncology}, pp. 363--370. New York: Plenum.

Wellner, J. A. and Zhang, Y. (2007) Two likelihood-based semiparametric estimation
methods for panel count data with covariates. \emph{Annals of Statistics}, \bold{35}(5), 2106--2142.

Lu, M., Zhang, Y. and Huang, J. (2009) Semiparametric estimation methods for panel
count data using monotone B-Splines. \emph{Journal of the American Statistical Association}
\bold{104}(487), 1060--1070.
}
\seealso{
\code{\link{blaTum}}
}
