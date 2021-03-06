% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parsers.R
\name{writeData}
\alias{writeData}
\title{writes out all data from the `platexpress` object}
\usage{
writeData(data, file, yids, xid, dtype = "data", verb = 1, ...)
}
\arguments{
\item{data}{\code{platexpress} data object}

\item{file}{either a character string naming a file or a connection open
for writing.  ‘""’ indicates output to the console.}

\item{yids}{IDs of the data to be written; if missing, all data will
written}

\item{xid}{ID of a data-set in the input data that can be used as x-axis
instead of the default Time vector}

\item{dtype}{type of the data to be plotted, default is the main 'data', but
e.g. 'orig' allows to plot the original data w/o processing (e.g.
interpolation, blank correction, etc.)}

\item{verb}{print messages if \code{>1}}

\item{...}{arguments passed to write.csv}
}
\description{
calls \code{\link{write.csv}} for each data type and writes
them to separate files names `<file>_<dataID>.csv`.
}
