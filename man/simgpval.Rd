\name{simgpval}
\alias{simgpval}
\title{Simulates Gaussian P-values}
\usage{
simgpval(y,x,i,nsim, plt=TRUE)
}
\arguments{
\item{y}{Dependent variable}

\item{x}{Covariates}

\item{i}{The chosen covariate}

\item{nsim}{Number of simulations}

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
a<-flag(snspt,3253,1,12)
simgpval(a[[1]],a[[2]],7,10,plt=FALSE)
}

