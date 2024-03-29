% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fgr1st.R
\name{fgr1st}
\alias{fgr1st}
\title{Calculates a dependence graph using Gaussian stepwise selection}
\usage{
fgr1st(x,p0=0.01,ind=0,kmn=0,kmx=0,mx=21,nedge=10^5,inr=T,xinr=F,qq=-1)
}
\arguments{
\item{x}{The matrix of covariates}

\item{p0}{Cut-off P-value}

\item{ind}{Restricts the dependent nodes to this subset}

\item{kmn}{The minimum number selected variables for each node irrespective of cut-off P-value}

\item{kmx}{The maximum number selected variables for each node irrespective of cut-off P-value}

\item{mx}{Maximum number of selected covariates for each node for all subset search} 

\item{nedge}{Maximum number of edges}

\item{inr}{Logical, if TRUE include an intercept}

\item{xinr}{Logical, if TRUE intercept already included}

\item{qq}{The number of covariates to choose from. If qq=-1 the number of covariates of x is used}


}
\value{
ned Number of edges

edg List of edges together with P-values for each  edge and proportional reduction of sum of squared residuals.
}
\description{
Calculates an independence graph using Gaussian stepwise selection
}
\examples{
data(boston)
a<-fgr1st(boston[,1:13],ind=3:6)
}
