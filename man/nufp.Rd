\name{nufp}
\docType{data}
\alias{nufp}
\title{nufp}
\description{This data set gives the results of simulations on the average number of false positives when the covariates are independent Gaussian noise. It is a 280 x 12 matrix.The first column is the number of covariates k, the second is the sample size n, columns 3:12 give the average number of false positives for the values of the order statistic nu for nu:1:10. The first 70 rows are for p0=0.01, the rows 71:140 p0=0.05, the rows 141:210 for p0=0.01/k , the rows 211:280 for p0=0.05/k
}

\usage{nufp}
\format{A 280 x 12 matrix.}
\source{Author's simulations}
\references{
"Covariate Selection Based on a Model-free Approach to Linear Regression with Exact Probabilities", 2020, <arXiv:1906.01990>
}
\keyword{datasets}