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
y The ith covariate of x without a lag, the dependent variable.

xl The covariates with lags from 1 :lag starting with the first covariate.
}
\description{
Calculation of lagged covariates
}
\examples{
data(vardata)
lvardata<-flag(vardata,256,1,182)
a<-f1st(lvardata[[1]],lvardata[[2]])
}
