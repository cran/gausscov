\name{f3sti}
\alias{f3sti}
\title{Selection of covariates with given excluded covariates}
\usage{
f3sti(y,x,covch,ind,m,p0=0.01,kmn=0,kmx=0, kex=0,mx=21,sub=T,inr=F,xinr=F,qq=-1,kexmx=100)
}
\arguments{
\item{y}{Dependent variable}

\item{x}{Covariates}

\item{covch}{Sum of squared residuals and selected covariates}

\item{ind}{The excluded covariates}

\item{m}{Number of iterations}



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
ind1 The excluded covariates

covch The sum of squared residuals and the selected covariates ordered in increasing size of sum of squared residuals

}
\description{
Selection of covariates with given excluded covariates
}
\examples{
data(leukemia)
covch=c(2.023725,1182,1219,2888,0)
covch<-matrix(covch,nrow=1,ncol=5)
ind<-c(1182,1219,2888)
ind<-matrix(ind,nrow=3,ncol=1)
m<-1
a<-f3sti(leukemia[[1]],leukemia[[2]],covch,ind,m,kexmx=5)
}
