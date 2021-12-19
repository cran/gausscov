#'  Calculation of dependence graph using repeated Gaussian stepwise selection.
#'
#' @param x Matrix of covariates.
#' @param p0  The cut-off P-value.
#' @param ind Restricts the dependent nodes to this subset
#' @param nu The order statistic of Gaussian covariates used for comparison.
#' @param kmn The minimum number of selected covariates for each node irrespective of cut-off P-value.
#' @param kmx The maximum number of selected covariates for each node irrespective of cut-off P-value.
#' @param nedge The maximum number of edges.
#' @param inr Logical, if TRUE  to include intercept.
#' @param xinr Logical, if TRUE intercept already included.
#' @return ned Number of edges
#' @return edg List of edges giving nodes (covariates), the approximations for each node, the covariates in the approximation and the corresponding P-values.
#' data(redwine)
#' a<-fgr2st(redwine[,1:11],ind=4:8) 
fgr2st<-function(x,p0=0.01,ind=0,nu=1,kmn=0,kmx=0,nedge=10^5,inr=T,xinr=F){
	n<-length(x[,1])
	k<-length(x)/n
	ind<-matrix(ind,nrow=1)
	if(ind[1]==0){ind<-1:k}
	li<-length(ind)
	p0<-p0/li
	if(inr){tmpx<-double(n)+1
		x<-cbind(x,tmpx)
		k<-k+1
		inr<-F
		xinr<-T
	}
	x<-matrix(x,nrow=n)
	xx<-x
	if(kmx==0){kmx<-min(n,k)}
	tmp<-.Fortran(
		"graphstst",
		as.double(xx),
		as.double(x),
		as.integer(n),
		as.integer(k),
		double(n),
		double(n),
		double(n),
		integer(k),
		as.double(p0),
		as.integer(kmx),	
		double((k+1)*2),
		integer(nedge*3),
		integer(1),
		integer(k),
		as.logical(xinr),
		as.double(nu),
		double(k),
		as.integer(nedge),
		double(k),
		double(nedge),
		as.integer(kmn),
		as.integer(li),
		as.integer(ind)
		)
	ned<-tmp[[13]]
	if(ned>0){
		edg<-tmp[[12]]
		edg<-matrix(edg,ncol=3)
		edg<-cbind(edg,tmp[[20]])
		edg<-edg[1:ned,]
	}
	else{
		edg<-NaN
	}
	list(ned,edg)
}
