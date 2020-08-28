#'  Calculation of dependency graph using Gaussian stepwise selection.
#'
#' @param x The variables
#' @param p0  Cut-off P-value 
#' @param km The maximum number of selected covariates for each node.
#' @param nu The order statistic of Gaussian covariates used for comparison.
#' @param nedge The maximum number of edges.
#' @param inr Logical if TRUE include an intercept.

#' @return ned Number of edges
#' @return edg List of edges
#' data(boston)
#' a<-fgr1st(boston[,1:13]) 
fgr1st<-function(x,p0=0.01,km=0,nu=1,nedge=10^6,inr=T,dr=F){
	n<-length(x[,1])
	k<-length(x)/n
	p0<-p0/k
	x<-matrix(x,nrow=n)
	xx<-x
	if(inr){tmpx<-double(n)+1
		x<-cbind(x,tmpx)
		x<-matrix(x,nrow=n)
		kk<-k+1
	}
	if(km==0){km=min(n,k)}
	km1<-km+1
	kexc<-integer(kk+1)
	tmp<-.Fortran(
		"graphst",
		as.double(xx),
		as.double(x),
		as.integer(n),
		as.integer(k),
		double(n),
		double(n),
		double(n),
		integer(kk+1),
		as.double(p0),
		as.integer(km),	
		double(km1*2),
		as.integer(km1),
		integer(k*km1*2),
		integer(1),
		as.integer(kexc),
		as.logical(inr),
		as.double(nu),	
		double(km1),
		as.integer(nedge),
		double(kk),
		as.integer(kk)
		)
	ned<-tmp[[14]]
	if(ned>0){
        edg<-tmp[[13]]
        edg<-matrix(edg,ncol=2)
        edg<-edg[1:ned,]
        if(!dr){
            tmp<-.Fortran(
                    "edge",
                    as.integer(edg),
                    as.integer(ned),
                    as.integer(ned),
                    integer(ned),
                    integer(1)
                )   
                ned<-tmp[[5]]
                edg<-tmp[[1]]
                edg<-matrix(edg,ncol=2)
        }
    }
    else{
        edg<-NaN
	}
	list(ned,edg)
}
