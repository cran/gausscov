\name{fvauto}
\alias{fvauto}
\title{Vector autoregressive approximation}
\usage{
fvauto(x,n,omx,p0=0.01,kmn=0,kmx=0,mx=21,kex=0,sub=T,inr=TRUE)
}
\arguments{
\item{x}{Variable}

\item{n}{Sample size}

\item{omx}{Maximum lag}

\item{p0}{The P-value cut-off}

\item{kmn}{Minimum number of included covariates irrespective of cut-off P-value}

\item{kmx}{Maxmum number of included covariates irrespective of cut-off P-value}

\item{mx}{The maximum number covariates for an all subset search}

\item{kex}{The excluded covariates}

\item{sub}{Logical, if TRUE best subset selected}

\item{inr}{Logical, if TRUE include intercept if not present}

}
\value{
res The selected lagged variables for each variable 
res2 The regression coefficients and P-values
res4 The residuals
}
\description{
Vector autoregressive approximation
}
\examples{
data(abcq)
 a<-fvauto(abcq,240,16)
}

