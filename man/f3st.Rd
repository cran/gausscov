\name{f3st}
\alias{f3st}
\title{Stepwise selection of covariates}
\usage{
f3st(y,x,m,mxa,p0=0.01,kmn=0,kmx=0,kex=0,mx=21,sub=T,inr=T,xinr=F,qq=-1,kexmx=100)}
\arguments{
\item{y}{Dependent variable}

\item{x}{Covariates}

\item{m}{The number of iterations}

\item{mxa}{The maximum number of approximations}

\item{p0}{The P-value cut-off}

\item{kmn}{The minimum number of included covariates irrespective of cut-off P-value}

\item{kmx}{The maximum number of included covariates irrespective of cut-off P-value.}

\item{kex}{The excluded covariates}

\item{mx}{The maximum number covariates for an all subset search}

\item{sub}{Logical if TRUE best subset selected}

\item{inr}{Logical if TRUE include intercept if not present}

\item{xinr}{Logical if TRUE intercept already present}

\item{qq}{The number of covariates to choose from. If qq=-1 the number of covariates of x is used.}

\item{kexmx}{The maximum number of covariates in an approximation.}

}
\value{
covch The standard deviation of the residuals and the selected covariates ordered in increasing size of the stamdard deviations.

lai  The number of rows of covch
}
\description{
Stepwise selection of covariates
}
\examples{
data(leukemia)
a<-f3st(leukemia[[1]],leukemia[[2]],m=2,mxa=7,kmn=5,sub=TRUE,kexmx=10)
}
