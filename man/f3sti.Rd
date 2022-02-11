\name{f3sti}
\alias{f3sti}
\title{Selection of covariates with given excluded covariates}
\usage{
f3sti(y,x,covch,ind,m,kexmx=100,p0=0.01,nu=1,kmn=0,kmx=0,

 mx=21,kex=0,sub=T,inr=T,xinr=F,qq=0)
}
\arguments{
\item{y}{Dependent variable}

\item{x}{Covariates}

\item{covch}{Sum of squared residuals and selected covariates}

\item{ind}{The excluded covariates}

\item{m}{Number of iterations}

\item{kexmx}{The maximum number of excluded covariates}

\item{p0}{The P-value cut-off}

\item{nu}{The order statistic of Gaussian covariates used for comparison}

\item{kmn}{The minimum number of included covariates irrespective of cut-off P-value}

\item{kmx}{The maximum number of included covariates irrespective of cut-off P-value.}

\item{mx}{The maximum number covariates for an all subset search}

\item{kex}{The excluded covariates}

\item{sub}{Logical if TRUE best subset selected}

\item{inr}{Logical if TRUE include intercept if not present}

\item{xinr}{Logical if TRUE intercept already present}

\item{qq}{The number of covariates to choose from. If qq=0 the number of covariates of x is used.}
}

\value{
ind1 The excluded covariates
covch The sum of squared residuals and the selected covariates ordered in increasing size of sum of squared residuals
}
\description{
Selection of covariates with given excluded covariates
}
\examples{
data(leukemia)
covch=c(2.023725,1182,1219,2888,0)
ind<-c(1182,1219,2888,0,0,0,0,0,0)
ind<-matrix(ind,ncol=3)
m<-3
a<-f3sti(leukemia[[1]],leukemia[[2]],covch,ind,m)
}
