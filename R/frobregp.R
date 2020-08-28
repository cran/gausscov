#' Robust regression using Huber's psi-function which calculates P-values
#'
#' @param y Dependent variable
#' @param x Covariates
#' @param cn  Tuning parameter for Huber's psi-function
#' @param sg  Scale
#' @param q  the numer of covariates available
#' @param ind The subset of covariates for which the results are required.
#' @param inr Logical TRUE to include intercept 
#' @param xinr Logical TRUE if x already has intercept.
#' @return ppi In order the subset ind, the regression coefficients, the P-values, the standard P-values.
#' @return res  Residuals
#' @return sg  Scale
#' @examples 
#' data(boston)
#' a<-frobregp(boston[,14],boston[,1:13])
frobregp<-function(y,x,cn=1,sg=0,q=-1,ind=0,inr=T,xinr=F){
	if(mad(y)==0){stop("MAD=0")}
	if(sg==0){sg<-mad(y)}
	tmpx<-cn*(1:1000)/1000
	cpp<-sum(tmpx^2*dnorm(tmpx))*cn/1000+cn**2*(1-pnorm(cn))
	cpp<-2*cpp
	n<-length(y)
	k<-length(x)/n
	x<-matrix(x,nrow=n)
	y<-matrix(y,ncol=1)
	if((!xinr)&inr){tmpx<-double(n)+1
		x<-cbind(x,tmpx)
		x<-matrix(x,nrow=n)
		xinr<-T
		k<-length(x)/n
	}
	q<-max(q,k)
	ind<-matrix(ind,ncol=1)
	if(ind[1]>0){
		if(!xinr){x<-x[,ind]}
		if(xinr){ind<-c(ind,k)
			x<-x[,ind]
		}
	}
	else{ind<-1:k}
	kk<-length(x)/n
	tmp<-.Fortran(	
		"robregp",
		as.double(y),
		as.double(x),
		double(n),
		double(n*kk),
		double(n*kk),
		double(kk^2),
		as.integer(n),
		as.integer(kk),
		double(kk),
		double(kk),
		double(kk),	
		double(n),	
		double(kk),
		as.double(cn),
		as.double(sg),
		double(3),
		as.double(cpp),
		integer(kk),
		double(2*kk),
		as.integer(q),
		as.logical(xinr)
		)
	beta<-tmp[[11]]
	beta<-matrix(beta,ncol=1)
	pp<-tmp[[19]]
	pp<-matrix(pp,ncol=2)
	res<-tmp[[12]]	
	sg<-tmp[[15]]
	ind<-matrix(ind,ncol=1)
	if(xinr){ind[kk]<-0}
	ppi<-cbind(ind,beta,pp)
	ppi<-matrix(ppi,ncol=4)
	list(ppi,res,sg)
}