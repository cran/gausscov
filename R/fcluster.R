#' Determine the disjoint connected components of an undirected dependency graph
#'
#' @param edg The edges of the graph.
#' @param q The number of covariates of the initial data from which the graph was obtained.
#' @return ncomp The number of components.
#' @return szcomp The sizes of the components.
#' @ return comp The covriates forming the components with alternating sign.
fcluster<-function(edg,q){
	nedg<-length(edg[,1])
	nodes<-unique(c(edg[,1],edg[,2]))
	nnodes<-length(nodes)
	tmp1<-unique(edg[,1])
	endnodes<-integer(q)+1
	endnodes[tmp1]<-0
	tmp<-.Fortran(
		"cluster",
		as.integer(edg),
		as.integer(endnodes),
		integer(q),
		integer(nnodes),
		integer(nnodes),
		integer(nnodes),
		as.integer(nedg),
		as.integer(nnodes),
		as.integer(q),
		integer(nnodes),
		integer(1)
		)
	ncomp<-tmp[[11]]
	szcomp<-tmp[[10]][1:ncomp]
	comp<-tmp[[4]]
	list(ncomp,szcomp,comp)
}
	
