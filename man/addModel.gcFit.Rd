% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interfaces.R
\name{addModel.gcFit}
\alias{addModel.gcFit}
\title{Add results from  \code{\link{grofit_plate}}
to \code{platexpress} object.}
\usage{
\method{addModel}{gcFit}(fit, data, ID = "model", ...)
}
\arguments{
\item{fit}{a \code{\link[grofit:gcFit]{gcFit}} object}

\item{data}{\code{platexpress} data object}

\item{ID}{data ID for the new object}

\item{...}{arguments passed to \code{\link{addData}}}
}
\description{
Calls the \code{\link[stats:predict]{predict}} method
for the growth curves fits returned
by grofit's \code{\link[grofit:gcFit]{gcFit}} or via \code{platexpress}'s
\code{\link{gcFit.2}} or
\code{\link{grofit_plate}}, and adds it to the
\code{platexpress} data object.
}
