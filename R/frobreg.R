#' Robust regression using Huber's psi-function or Hampel's redescending psi-function without P-values
#'
#' @param y Dependent variable
#' @param x Covariates
#' @param cn  Tuning parameter for Huber's psi-function
#' @param cnr The constants for Hampel's three part redescending psi function
#' @param sg  Scale
#' @param scale Logical, if TRUE calculate sg simultaneously, otherwise keeps initial sg
#' @param inr Logical if TRUE to include intercept 
#' @param xinr Logical if TRUE intercept already included
#' @param red  Logical If true Hampel's three part redescending psi function
#' @return beta Regression coefficients
#' @return res  Residuals
#' @return sg  Scale
#' @return rho Sums of rho, psi and psi1 functions.
#' @examples 
#' data(boston)
#' a<-frrg(boston[,14],boston[,1:13])
frrg<-function(y,x,cn=1,cnr=c(2,4,8),sg=0,scale=T,inr=T,xinr=F,red=F){
	if(mad(y)==0){stop("MAD=0")}
	cnn<-cn
	if(red){cnn<-cnr[1]}
	tmpx<-cnn*(1:1000)/1000
	cpp<-sum(tmpx^2*dnorm(tmpx))*cnn/1000+cnn**2*(1-pnorm(cnn))
	cpp<-2*cpp
	n<-length(y)
	x<-matrix(x,nrow=n)
	y<-matrix(y,ncol=1)
	if((!xinr)&inr){tmpx<-double(n)+1
		x<-cbind(x,tmpx)
		x<-matrix(x,nrow=n)
		xinr<-T
	}
	if(xinr){mny<-median(y)
		y<-y-mny
	}
	if(sg==0){sg<-median(abs(y))}
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
		as.logical(scale),
		as.logical(red),
		as.double(cnr)
		)
	beta<-tmp[[10]]
	res<-tmp[[11]]	
	sg<-tmp[[14]]
	rho<-tmp[[15]]
	if(xinr){beta[k]<-beta[k]+mny}
	list(beta,res,sg,rho)
}

