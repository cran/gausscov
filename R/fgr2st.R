#'  Calculation of dependency graph using repeated Gaussian stepwise procedure.
#'
#' @param x Matrix of covariates
#' @param p0  The P-value cut-off
#' @param km The maximum number of included covariates for each covariate.
#' @param nu The order statistic of Gaussian covariates used for comparison.
#' @param nedge The maximum number of edges.
#' @param inr Logical, if TRUE  to include intercept.
#' @param dr  Logical, if TRUE (a,b) and (b,a), a not equal b,  are different edges
#' @return ned Number of edges
#' @return edg List of edges
#' data(redwine)
#' a<-fgr2st(redwine[,1:11]) 
fgr2st<-function(x,p0=0.01,km=0,nu=1,nedge=10^6,inr=T,dr=F){
	n<-length(x[,1])
	k<-length(x)/n
	p0<-p0/k
	x<-matrix(x,nrow=n)
	xx<-x
	kk<-k
	if(inr){tmpx<-double(n)+1
		dim(tmpx)<-c(n,1)
		x<-cbind(x,tmpx)
		kk<-k+1
	}
	if(km>0){
		if(inr){km<-km+1}
	}
	if(km==0){km=min(n,k)}
	km1<-km+1
	kexc<-integer(kk)+1
	km1<-km+1
	nee<-k*km1*3
	grph<-integer(nee)
	grph<-matrix(grph,ncol=3)
	tmp<-.Fortran(
		"graphstst",
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
		as.integer(grph),
		integer(1),
		as.integer(kexc),
		as.integer(nedge),
		as.logical(inr),
		as.double(nu),
		double(km1),
		double(kk),
		as.integer(kk)
		)
	ned<-tmp[[14]]
	if(ned>0){
		edg<-tmp[[13]]
		edg<-matrix(edg,ncol=3)
		edg<-edg[1:ned,]
		edg<-matrix(edg,ncol=3)
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
