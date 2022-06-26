\name{fcluster}
\alias{fcluster}
\title{Disjoint components of an undirected dependency graph}
\usage{
fcluster(edg,q)
}
\arguments{
\item{edg}{the edges of the graph}

\item{q}{The number of covariates used to construct the graph}

}
\value{
ncomp The number of components.

szcomp The sizes of the components

comp The covriates forming the components with alternating sign.
}
\description{
Determine the disjoint connected components of an undirected dependency graph
}
\examples{
data(boston)
bostint<-fgeninter(boston[,1:13],2)[[1]]
a<-f1st(boston[,14],bostint,kmn=10,sub=TRUE)
}

