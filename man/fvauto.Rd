\name{fvauto}
\alias{fvauto}
\title{Vector autoregressive approximation}
\usage{
fvauto(x,n,omx,p0=0.01)
}
\arguments{
\item{x}{Variable}

\item{n}{Sample size}

\item{omx}{Maximum lag}

\item{p0}{The P-value cut-off}

}
\value{
res The selected lagged variables for each variable 
}
\description{
Vector autoregressive approximation
}
\examples{
data(abcq)
 a<-fvauto(abcq,240,10)
}

