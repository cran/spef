\name{panelReg.object}
\alias{panelReg.object}
\title{Semiparametric Panel Count Regression Object}
\description{
  This S3 class of objects is returned by the \code{panelReg} class of
  functions to represent a fitted semiparametric panel count regression
  model. Objects of this class have methods for the functions
  \code{print} and \code{plot}.
}
\arguments{
  \item{beta}{a vector, estimated coefficients of the linear predictor.}
  \item{baseline}{a step function, estimated cumulative baseline mean.}
  \item{lambda}{a vector, estimated baseline mean for each interval}
  \item{convergence}{an integer code, \code{0} indicates successful
    convergence, error code \code{1} indicates that the iteration limit
    \code{maxIter} had been reached.}
  \item{iter}{an integer value, number of interations when stopped.}
  \item{betaSE}{a vector, estimated standard errors of beta.}
  \item{betaVar}{a matrix, estimated covariance matrix of beta.}
  \item{baselineSE}{a step function, estimated standard error of
    cumulative baseline mean.}
}
\seealso{
  \code{\link{panelReg}}
}
\keyword{survival}
