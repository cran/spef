\name{blaTum}
\docType{data}
\alias{blaTum}
\title{Bladder Tumors Cancer Recurrences}
\usage{
blaTum
}
\description{A data frame contains data on recurrences of bladder
  cancer, used by many people to demonstrate methodology for recurrent
  event modelling. It's used in this package to illustrate methodology
  for panel count modeling.
}
\format{
  \tabular{ll}{
    id:\tab patient id\cr
    treatment:\tab placebo = 0, thiotepa = 1\cr
    size:\tab size (cm) of largest initial tumour\cr
    num:\tab initial number of tumours (8=8 or more)\cr
    time:\tab observation time\cr
    count:\tab number of new tumors since last observation time\cr
  }
}
\source{
  Byar, D. P. (1980). The Veterans administration study of
  chemoprophylaxis for recurrent stage I bladder tumors: comparisons
  of placebo, pyridoxine, and topical thiotepa. \emph{Bladder Tumors and
  Other Topics in Urological Oncology}, pp. 363--370. New York: Plenum.

  Huang, C. Y., Wang, M. C. and Zhang, Y. (2006). Analysing panel count data
  with informative observation times. \emph{Biometrika}, \bold{93}(4):763--776.
}
\examples{
# Plot bladder tumor data
ggplot(blaTum, aes(time, id)) + geom_tile(aes(fill=count)) +
   facet_grid(treatment ~ ., scales="free_y", )

# Fit a semiparametric regression model
formula <- PanelSurv(id, time, count) ~ treatment + num + size
fit <- panelReg(formula, data=blaTum, method="AEE", se="Sandwich")
fit
}
\keyword{survival}
\keyword{panel count}
