\name{f3st}
\alias{f3st}
\title{Stepwise selection of covariates}
\usage{
f3st(y,x,m,kexmx=100,p0=0.01,nu=1,kmn=0,kmx=0,mx=21,lm=1000,kex=0,sub=T,inr=T,xinr=F,qq=0)
}
\arguments{
\item{y}{Dependent variable}

\item{x}{Covariates}

\item{m}{The number of iterations}

\item{kexmx}{The maximum number of covariates in an approximation}

\item{p0}{The P-value cut-off}

\item{nu}{The order statistic of Gaussian covariates used for comparison}

\item{kmn}{The minimum number of included covariates irrespective of cut-off P-value}

\item{kmx}{The maximum number of included covariates irrespective of cut-off P-value.}

\item{mx}{The maximum number covariates for an all subset search}

\item{lm}{The maximum number of approximations.}

\item{kex}{The excluded covariates}

\item{sub}{Logical if TRUE best subset selected}

\item{inr}{Logical if TRUE include intercept if not present}

\item{xinr}{Logical if TRUE intercept already present}

\item{qq}{The number of covariates to choose from. If qq=0 the number of covariates of x is used.}
}
\value{
covch The sum of squared residuals and the selected covariates ordered in increasing size of sum of squared residuals.

lai  The number of rows of covch
}
\description{
Stepwise selection of covariates
}
\examples{
data(leukemia)
a<-f3st(leukemia[[1]],leukemia[[2]],m=2,kexmx=5,kmn=5,sub=TRUE)
}
