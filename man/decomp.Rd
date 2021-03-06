% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fdecomp.R
\name{decomp}
\alias{decomp}
\title{decompose a given interaction ic into its component parts}
\usage{
decomp(ic, k, ord)
}
\arguments{
\item{ic}{The number of the interaction}

\item{k}{The number of covariates of x without intercept} 

\item{ord}{The order of the interactions}
}
\value{
decom The component parts of the interaction.
}
\description{
decompose a given interaction ic into its component parts
}
\examples{
a<-decomp(7783,14,8)
}
