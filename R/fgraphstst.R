#' Calculates an independence graph for a set of variables using  repeated stepwise selection 
#'
#' @param x  The variables
#' @param alpha Cut-off p-value
#' @param nu Order statistic
#' @param kmax Maximum number selected variables for each node
#' @param nedge Maximum number of edges
#' @param intercept If true intercept included
#' @param chkintercept If true intercept included depending on p-value
#' @return ned Number of edges
#' @return edg The edges for each node in the graph
#' @examples 
#' colgrph<-fgraphstst(colon.x,0.05)
#' colgrph[[1]]
#' colgrph[[2]][1:20,]
fgraphstst<-function(x,alpha,nu=1,kmax=0,nedge=10^6,intercept=TRUE,chkintercept=FALSE){
	k<-length(x[1,])
	alpha<-alpha/k
	n<-length(x[,1])
	x<-as.matrix(x)
	if(kmax==0){kmax<-n}
	if(intercept){
		tmpx<-double(n)+1
		dim(tmpx)<-c(n,1)
		x<-cbind(tmpx,x)
		k<-k+1
		kmax<-kmax+1
	}
	kmax1<-kmax+1
	kexc<-integer(k)
	tmp<-.Fortran(
		"graphstst",
		as.double(x),
		as.double(x),
		as.integer(n),
		as.integer(k),
		double(n),
		double(n),
		double(n),
		integer(k),
		as.double(alpha),
		as.integer(kmax),	
		double(kmax1*2),
		as.integer(kmax1),
		integer(nedge*2),
		integer(1),
		as.integer(kexc),
		as.integer(nedge),
		as.logical(intercept),
		as.double(nu),
		double(kmax1),
		integer(k+1),
		integer(1),
		as.logical(chkintercept),
		integer(kmax1)
		)
	if(tmp[[16]]>nedge) stop("graph too large")
	edg<-tmp[[13]]
	if(max(edg)>0){
		dim(edg)<-c(nedge,2)
		ne<-tmp[[14]]
		edg<-edg[1:ne,]
		kmax<-1000
		tmp<-.Fortran(
			"edge",
			as.integer(edg),
			as.integer(ne),
			as.integer(kmax),
			integer(kmax),
			integer(1)
		)
		edg<-tmp[[1]]
		dim(edg)<-c(ne,2)
		ned<-tmp[[5]]
		edg<-edg[1:ned,]
	}
	else{
		ned<-0
		edg<-0
	}
	list(ned,edg)
}
