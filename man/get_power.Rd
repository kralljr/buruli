% Generated by roxygen2 (4.0.1): do not edit by hand
\name{get_power}
\alias{get_power}
\title{Seasonal power}
\usage{
get_power(N, a, w, P, t, sd = 1, sig = 0.05, max1 = 0)
}
\arguments{
\item{N}{Number of simulations}

\item{a}{amplitude}

\item{w}{frequency}

\item{P}{phase}

\item{t}{time}

\item{sd}{standard deviation of random error}

\item{sig}{significance level (defaults to 0.05)}

\item{max1}{Maximum number of days delayed}
}
\description{
\code{get_power} Obtains power for seasonal associations using cosignor (Barnett and Dobson 2010)
}
\examples{
N <- 100
a <- 1
w <- 1/16
P <- 0
t <- seq(1, 365 * 2, by = 1)
get_power(100, a, w, P, t)
}

