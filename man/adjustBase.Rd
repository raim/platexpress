% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/platereader.R
\name{adjustBase}
\alias{adjustBase}
\title{Adjusts to Minimal Base}
\usage{
adjustBase(data, base = 0, yids, add.fraction, xlim, each = FALSE, verb = TRUE)
}
\arguments{
\item{data}{\code{\link{platexpress}} data, see \code{\link{readPlateData}}}

\item{base}{the new minimum for the data, default is 0, but it could
e.g. be the OD used for inoculation}

\item{yids}{vector of ID strings for which base correction should be
executed}

\item{add.fraction}{a fraction of the whole data range, added to base}

\item{xlim}{min and max row number of the data to be adjusted}

\item{each}{add base for each well separately!}

\item{verb}{print messages if true}
}
\value{
Returns `data' where all data sets (or only those specified in
option \code{yids}) where raised to a minimum level in
}
\description{
The function raises all data to a "base" level, default 0, to avoid
negative values that sometimes occur after blank correction
in \code{\link{correctBlanks}}. The function can be optionally
called directly in \code{\link{correctBlanks}} by option \code{base}.
}
\details{
Adjusts data to a new mininum, this is useful for adjustment
of negative values after blank corrections.
}
\seealso{
\code{\link{correctBlanks}}
}
\author{
Rainer Machne \email{raim@tbi.univie.ac.at}
}
