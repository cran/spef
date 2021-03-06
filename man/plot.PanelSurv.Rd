% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PanelSurv.R
\name{plot.PanelSurv}
\alias{plot.PanelSurv}
\title{Produce Tile Plot}
\usage{
\method{plot}{PanelSurv}(x, heat = FALSE, order = TRUE, ...)
}
\arguments{
\item{x}{an object of class \code{PanelSurv}.}

\item{heat}{an optional logical value indicating whether
a swimmer-plot-like tile plot will be produced.}

\item{order}{an optional logical value indicating whether the tile plot
will be sorted by the largest observationt time.}

\item{...}{future extension}
}
\value{
A \code{ggplot} object
}
\description{
Plot the tile plot from a \code{PanelSurv} object.
}
