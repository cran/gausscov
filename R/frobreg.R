#' Robust regression using Huber's psi-function
#'
#' @param y Dependent variable
#' @param x Covariates
#' @param cn  Tuning parameter for Huber's psi-function
#' @param sg  Scale
#' @param scale Logical, if TRUE calculate sg simultaneously, otherwise keeps initial sg
#' @param inr Logical if TRUE to include intercept 
#' @param xinr Logical if TRUE intercept already included
#' @return beta Regression coefficients
#' @return res  Residuals
#' @return sg  Scale
#' @examples 
#' data(boston)
#' a<-frobreg(boston[,14],boston[,1:13])
frobreg<-function(y,x,cn=1,sg=0,scale=T,inr=T,xinr=F){
	if(mad(y)==0){stop("MAD=0")}
	if(sg==0){sg<-mad(y)}
	tmpx<-cn*(1:1000)/1000
	cpp<-sum(tmpx^2*dnorm(tmpx))*cn/1000+cn**2*(1-pnorm(cn))
	cpp<-2*cpp
	n<-length(y)
	x<-matrix(x,nrow=n)
	y<-matrix(y,ncol=1)
	if((!xinr)&inr){tmpx<-double(n)+1
		x<-cbind(x,tmpx)
		x<-matrix(x,nrow=n)
		xinr<-T
	}
	k<-length(x)/n
	tmp<-.Fortran(	
		"robreg",	
		as.double(y),
		as.double(x),
		double(n),
		double(n*k),
		double(k^2),
		as.integer(n),
		as.integer(k),
		double(k),
		double(k),
		double(k),	
		double(n),	
		double(k),
		as.double(cn),
		as.double(sg),
		double(3),
		as.double(cpp),
		as.logical(scale)
		)
	beta<-tmp[[10]]
	res<-tmp[[11]]	
	sg<-tmp[[14]]
	list(beta,res,sg)
}

