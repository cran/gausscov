\name{flag}
\alias{flag}
\title{Calculation of lagged covariates}
\usage{
flag(x,n,lag)
}
\arguments{
\item{x}{The covariates}

\item{n}{Sample size}

\item{lag}{Themaximum lag}
}
\value{
y The first covariate without a lag, the dependent covariate.

xl The covaraites with lags from 1 :lag starting with the first covariate.
}
\description{
Calculation of lagged covariates
}
\examples{
data(abcq)
abcql<-flag(abcq,240,16)
a<-f1st(abcql[[1]],abcql[[2]])
}
