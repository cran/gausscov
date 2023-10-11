\name{decode}
\alias{decode}
\title{decodes subsets}
\usage{
decode(ns,k)
}
\arguments{
\item{ns}{The number of the subset}

\item{k}{The number of covariates excluding the intercept}
}

\value{
ind the list of covariates

set  A binary vector giving the covariates
}
\description{
Decodes the number of a subset selected by fasb.R to give the  covariates
}
\examples{
decode(19,8)
}

