#' Robust regression using Huber's psi-function or Hampel's three part redescending psi function which calculates P-values
#'
#' @param y Dependent variable.
#' @param x Covariates.
#' @param cn  Tuning parameter for Huber's psi-function.
#' @param cnr The constants for Hampel's three part redescending psi function.
#' @param sg  Scale.
#' @param q   The numer of covariates available.
#' @param ind  The subset of covariates for which the regression is required.
#' @param scale Logical TRUE to recalculate scale sg.
#' @param ind The subset of covariates for which the results are required.
#' @param inr Logical TRUE to include intercept.
#' @param xinr Logical TRUE if x already has intercept.
#' @param red  Logical It true Hampel's three part redescending psi function.
#' @return ppi In order the subset ind, the regression coefficients, the P-values, the standard P-values.
#' @return res  Residuals.
#' @return sg  Scale.
#' @return rho Sums of rho, psi and psi1 functions.
#' @examples 
#' data(boston)
#' a<-frrgp(boston[,14],boston[,1:13])
frrgp<-function(y,x,cn=1,cnr=c(1,2,4),sg=0,q=-1,ind=0,scale=T,inr=T,xinr=F,red=F){
	if(mad(y)==0){stop("MAD=0")}
	cnn<-cn
	if(red){cnn<-cnr[1]}
	tmpx<-cnn*(1:1000)/1000
	cpp<-sum(tmpx^2*dnorm(tmpx))*cnn/1000+cnn**2*(1-pnorm(cnn))
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
	ind<-matrix(ind,ncol=1)
	if(ind[1]>0){
		if(!xinr){x<-x[,ind]}
		if(xinr){ind<-c(ind,k)
			x<-x[,ind]
		}
	}
	else{ind<-1:k}
	kk<-length(x)/n
	if(xinr){mny<-median(y)
		y<-y-mny
	}
	if(sg==0){sg<-median(abs(y))}
	if(q==-1){q<-kk}
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
		as.logical(xinr),
		double(n),
		as.logical(scale),
		as.logical(red),
		as.double(cnr)
		)
	beta<-tmp[[11]]
	beta<-matrix(beta,ncol=1)
	if(xinr){beta[kk]<-beta[kk]+mny}
	pp<-tmp[[19]]
	pp<-matrix(pp,ncol=2)
	res<-tmp[[22]]	
	sg<-tmp[[15]]
	ind<-matrix(ind,ncol=1)
	if(xinr){ind[kk]<-0}
	ppi<-cbind(ind,beta,pp)
	ppi<-matrix(ppi,ncol=4)
	rho<-tmp[[16]]
	list(ppi,res,sg,rho)
}
