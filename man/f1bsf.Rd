\name{f1bsf}
\alias{f1bsf}
\title{Stepwise selection of interval covariates in non-paametric regression}
\usage{
f1bsf(y,k,lam,pr=0.5,mm=20,p0=0.01)
}
\arguments{
\item{y}{Dependent variable}

\item{k}{Length odf smallest interval}

\item{lam}{Factor for increasing  size of intervals}

\item{pr}{Proportional incerease in size of sample to reduce edge effects}

\item{mm}{Parameter of fgentrig for the number of trigonometric functions}

\item{p0}{Gaussian P-value threshold}
}

\value{
pv ff The approximation

res The residuals
}
\description{
Stepwise selection of interval covariates in non-paametric regression
}
\examples{
data(vardata)
a<-f1bsf(vardata[,70],4,1.1,pr=0,mm=10)
}
