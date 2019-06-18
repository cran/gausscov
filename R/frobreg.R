#' Robust regression using Huber's psi-function
#'
#' @param y Dependent variable
#' @param x Covariates
#' @param cn  Tuning parameter for Huber's psi-function
#' @param cpp Fisher consistency parameter
#' @param sig  Scale
#' @param intercept Logical to include intercept 
#' @return beta Regression coefficients
#' @return res  Residuals
#' @return sig  Scale
#' @examples 
#' frobreg(boston[,14],boston[,1:13],1)
frobreg<-function(y,x,cn,cpp=0,sig=0,intercept=TRUE){			
	n<-length(y)
	xx<-as.matrix(x)
	res<-y
	yy<-y
	if(cpp==0){
		tmpx<-cn*(1:1000)/1000
		cpp<-sum(tmpx^2*dnorm(tmpx))*cn/1000+cn**2*(1-pnorm(cn))
		cpp<-2*cpp
	}
	if(intercept){
	 x0<-double(n)+1			# intercept		
	 xx<-cbind(x0,xx)
	}	
	k<-length(xx[1,])
	beta<-double(k)+0
	if(sig==0){sig<-mad(y)}
	tmp<-.Fortran(			# robust regression with all covariates 
		"robreg",		# tuning constant cn and scale sig
		as.double(y),
		as.double(xx),
		as.double(yy),
		as.double(xx),
		double(k^2),
		as.integer(n),
		as.integer(k),
		double(k),
		double(k),
		as.double(beta),	
		as.double(res),	
		double(k),
		as.double(cn),
		as.double(sig),
		double(3),
		as.double(cpp)
		)
	beta<-tmp[[10]]
	res<-tmp[[11]]	
	sig<-tmp[[14]]
	list(beta,res,sig)
}
#
