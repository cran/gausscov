\name{flag}
\alias{flag}
\title{Calculation of lagged covariates}
\usage{
flag(x,n,i,lag)
}
\arguments{
\item{x}{The covariates}
\item{n}{The sample size}
\item{i}{The dependent variable}
\item{lag}{The maximum lag}
}
\value{
y The ith covariate of x without a lag, the dependent varaible.

xl The covariates with lags from 1 :lag starting with the first covariate.
}
\description{
Calculation of lagged covariates
}
\examples{
data(abcq)
abcql<-flag(abcq,240,1,16)
a<-f1st(abcql[[1]],abcql[[2]])
}
