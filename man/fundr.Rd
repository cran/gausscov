\name{fundr}
\alias{fundr}
\title{Converts directed into an undirected graph}
\usage{
fundr(gr)
}
\arguments{
\item{gr}{A directed graph}
}
\value{
gr The undirected graph
}
\description{
Conversion of a directed graph into an undirected graph
}
\examples{
data(boston)
grb<-fgr1st(boston[,1:13])
grbu<-fundr(grb[[2]][,1:2])
}
	