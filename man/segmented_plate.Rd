% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interfaces.R
\name{segmented_plate}
\alias{segmented_plate}
\title{Call \code{\link[segmented]{segmented}} for selected wells.}
\usage{
segmented_plate(
  data,
  yid = "OD",
  wells,
  log = TRUE,
  xid,
  man = 1,
  maxtry = 5,
  control = segmented::seg.control(),
  psis,
  psi,
  npsi = 5,
  plot = FALSE,
  verb = 0,
  ...
)
}
\arguments{
\item{data}{\code{platexpress} data object}

\item{yid}{ID of the \code{platexpress} data to use}

\item{wells}{subset and plot order of wells}

\item{log}{use ln of the data}

\item{xid}{x-axis data ID in \code{platexpress} data}

\item{man}{apply a moving average \code{\link{ma}} with \code{n=man}}

\item{maxtry}{maximum number of attempts to run
\code{\link[segmented]{segmented}}}

\item{control}{a list of parameters for controlling the
\code{\link[segmented]{segmented}} fitting process.
See the documentation for \code{\link[segmented:seg.control]{seg.control}}
for details.}

\item{psis}{named list of breakpoints for wells, generated from
argument \code{psi} if missing}

\item{psi}{vector of breakpoints (x-values) to be passed to
\code{\link[segmented]{segmented}}, generated as equal spaced breakpoints
from argument \code{npsi} if missing}

\item{npsi}{number of equally spaced breakpoints, used if both
\code{psis} and \code{psi} are missing}

\item{plot}{plot each fit}

\item{verb}{progress messages}

\item{...}{arguments passed to \code{\link[segmented]{segmented}}}
}
\description{
The package \code{\link[segmented]{segmented}} splits curves
into linear segments, given a list of pre-defined breakpoints.
If the passed data is a biomass measure
(eg. OD), and option \code{log=TRUE} the slopes correspond to local
growth rate (per time unit). Missing (NA) or non-finite
y-values will be removed. \code{\link[segmented]{segmented}} has
a random initialization and sometimes fails. The function will
attempt \code{maxtry} of calls before giving up.
}
