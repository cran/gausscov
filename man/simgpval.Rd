\name{simgpval}
\alias{simgpval}
\title{Simulates Gaussian P-values}
\usage{
simgpval(y,x,i,nsim,qq=-1,plt=TRUE)
}
\arguments{
\item{y}{Dependent variable}

\item{x}{Covariates}

\item{i}{The chosen covariate}

\item{nsim}{The number of simulations}

\item{qq}{The total number of covariates available}

\item{plt}{Logical, if TRUE the F P-values of the Gaussian covariates are ordered and plotted}

}
\value{
pvg P-value of x_i and relative frequency
}
\description{
Simulates Gaussian P-values
}
\examples{
data(snspt)
snspt<-matrix(snspt,nrow=3253,ncol=1)
a<-flag(snspt,3253,1,12,inr=FALSE)
simgpval(a[[1]],a[[2]],7,10,plt=FALSE)
}

