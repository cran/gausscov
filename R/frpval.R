#' Robust regression using Huber's psi-function or Hampel's three part redescending psi function providing P-values
#'
#' @param y Dependent variable.       
#' @param x Covariates.
#' @param ind The subset of covariates for which the results are required
#' @param cn  Tuning parameter for Huber's psi-function.
#' @param cnr The constants for Hampel's three part redescending psi function.
#' @param sg  Scale. If 0 the MAD is used.
#' @param q   The number of covariates available. If q=-1 the covariates are used.
#' @param scale Logical If TRUE scale sg recalculated.
#' @param inr Logical TRUE to include intercept.
#' @param xinr Logical TRUE if x already has intercept.
#' @param red  Logical If TRUE Hampel's three part redescending psi function.
#' @return ppi In order the subset ind, the regression coefficients, the Gaussian P-values, the standard P-values.
#' @return res  Residuals.
#' @return sg  Scale.
#' @return rho Sums of rho, psi and psi1 functions.
#' @examples 
#' data(boston)
#' a<-frpval(boston[,14],boston[,1:13],1:6)
frpval<-function(y,x,ind,cn=1,cnr=c(1,3,5),sg=0,q=-1,scale=T,inr=T,xinr=F,red=F){
	if(mad(y)==0){stop("MAD=0")}
	cnn<-cn
	if(red){cnn<-cnr[3]
		tmpx<-cnn*(1:1000)/1000
		cpp<-0
		for(i in 1:1000){
			cpp<-cpp+fpsired(tmpx[i],cnr)^2*dnorm(tmpx[i])
		}
		cpp<-cpp*cnn/1000
	}
	else{tmpx<-cnn*(1:1000)/1000
		cpp<-sum(tmpx^2*dnorm(tmpx))*cnn/1000+cnn**2*(1-pnorm(cnn))
	}
	cpp<-2*cpp
	n<-length(y)
	kx<-length(x)/n
	x<-matrix(x,nrow=n)
	y<-matrix(y,ncol=1)
	li<-length(ind)
	if(q==-1){	q<-kx-li+1}
	if(!xinr){
		if(inr){
         			tmpx<-double(n)+1
         			x<-cbind(x,tmpx)
			ind<-c(ind,kx+1)
			kx<-kx+1
 			xinr<-TRUE
			
		}
     	}
	ind<-matrix(ind,nrow=1)
	if(xinr){
		if(max(ind)<kx){ind<-c(ind,kx)}
	}
	xx<-x[,ind]
	kk<-length(xx)/n
	if(sg==0){sg<-mad(y)}
	tmp<-.Fortran(	
		"robregp",
		as.double(y),
		as.double(xx),
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
