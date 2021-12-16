\name{flag}
\alias{flag}
\title{Calculation of lagged covariates}
\usage{
flag(x,n,lag)
}
\arguments{
\item{x}{The covariates}

\item{n}{The sample size}

\item{lag}{The maximum lag}
}
\value{
y The first covariate of x without a lag, the dependent covariate.

xl The covariates with lags from 1 :lag starting with the first covariate.
}
\description{
Calculation of lagged covariates
}
\examples{
data(abcq)
abcql<-flag(abcq,240,16)
a<-f1st(abcql[[1]],abcql[[2]])
}
