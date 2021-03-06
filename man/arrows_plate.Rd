% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{arrows_plate}
\alias{arrows_plate}
\title{Draw arrows of growth phases.}
\usage{
arrows_plate(
  segments,
  wells,
  col = "#000000",
  ncol = "#F81894",
  wcol,
  xlim,
  maxmu,
  minmu,
  lwd.max = 10,
  head.length = 0.02,
  add = TRUE,
  axis2 = !add,
  ...
)
}
\arguments{
\item{segments}{a list of segmentations for each well as returned
by \code{\link{dpseg_plate}} or \code{\link{segmented_plate}}}

\item{wells}{subset and plot order of wells}

\item{col}{default color for arrows from positive slope segments}

\item{ncol}{default color for arrows from negative slope segments}

\item{wcol}{well colors (named vector) for arrows and axis labels,
eg. from plate layout map}

\item{xlim}{correct left and right borders required for
results from \code{segmented_plate}}

\item{maxmu}{maximal slope (growth rate) for line width scaling}

\item{minmu}{minimal slope (growth rate) for line width scaling}

\item{lwd.max}{maximal line width}

\item{head.length}{arrow head length, argument \code{length} to
\code{\link{arrows}}}

\item{add}{add to existing plot}

\item{axis2}{add wells as y-axis tick labels}

\item{...}{arguments passed to \code{\link{arrows}}}
}
\description{
draw arrows of linear segments calculated by
\code{\link{dpseg_plate}} or \code{\link{segmented_plate}}
with local slope (growth rate) as scaling factor, optionally in
well colors.
}
