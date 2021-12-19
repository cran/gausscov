#'  Calculation of dependence graph using Gaussian stepwise selection.
#'
#' @param x Matrix of covariates.
#' @param p0  Cut-off P-value. 
#' @param ind Restricts the dependent nodes to this subset
#' @param nu The order statistic of Gaussian covariates used for comparison.
#' @param kmn The minimum number of selected covariates for each node irrespective of cut-off P-value.
#' @param kmx The maximum number of selected covariates for each node irrespective of cut-off P-value.
#' @param nedge The maximum number of edges.
#' @param inr Logical if TRUE include an intercept.
#' @param xinr Logical if TRUE intercept already included.
#' @return ned Number of edges
#' @return edg List of edges together with P-values for each  edge and proportional reduction of sum of squared residuals
#' data(boston)
#' a<-fgr1st(boston[,1:13],ind=3:6) 
fgr1st<-function(x,p0=0.01,ind=0,nu=1,kmn=0,kmx=0,nedge=10^5,inr=T,xinr=F){
	n<-length(x[,1])
	k<-length(x)/n
	ind<-matrix(ind,nrow=1)
	if(ind[1]==0){ind<-1:k}
	li<-length(ind)
	p0<-p0/li
	x<-matrix(x,nrow=n)
	xx<-x
	if(inr){tmpx<-double(n)+1
		x<-cbind(x,tmpx)
		k<-k+1
		xinr<-T
		inr<-F
	}
	if(kmx==0){kmx<-min(n,k)}
	tmp<-.Fortran(
		"graphst",
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
		integer(nedge*2),
		integer(1),
		integer(k),
		as.logical(xinr),
		as.double(nu),	
		double(k),
		as.integer(nedge),
		double(k),
		as.integer(kmn),
		as.integer(li),
		as.integer(ind)
		)
	ned<-tmp[[13]]
	if(ned>0){
        		edg<-tmp[[12]]
        		edg<-matrix(edg,ncol=2)
        		edg<-edg[1:ned,]
        		edg<-matrix(edg,ncol=2)
		p<-double(2*ned)
		p<-matrix(p,ncol=2)
		for(i in 1:k){
			tmpi<-(1:ned)[edg[,1]==i]
			ind<-edg[tmpi,2]
			li<-length(ind)
			if(li>0){
				qq<-k-li
				if(xinr){qq<-qq-1}
				a<-fpval(x[,i],x,ind,q=qq,inr=F,xinr=xinr)
				ka<-length(a[[1]][,1])
				ika<-(1:ka)[a[[1]][,1]>0]
				p[tmpi,1]<-a[[1]][ika,3]
				p[tmpi,2]<-a[[1]][ika,5]
			}
		}
		edg<-cbind(edg,p)
    	}
    	else{
       		edg<-NaN
	}
	list(ned,edg)
}
