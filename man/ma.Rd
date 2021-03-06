% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/platereader.R
\name{ma}
\alias{ma}
\title{Moving Average}
\usage{
ma(x, n = 5, circular = FALSE)
}
\arguments{
\item{x}{data vector along which a moving average will be calculated}

\item{n}{moving average window size}

\item{circular}{logical see help of function \code{\link[stats]{filter}}}
}
\description{
calculate a moving average using \code{\link[stats]{filter}}
}
