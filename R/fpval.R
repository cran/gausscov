#' Calculates the regression coefficients, the P-values and the standard P-values for the chosen subset ind
#
#' @param y The dependent variable
#' @param x  The covariates
#' @param ind The indices of the subset whose P-values are required
#' @param q   The total  number of covariates from which ind was selected
#' @param inr Logical If TRUE intercept to be included
#' @param xinr Logical If TRUE intercept already included, overides inr
#' @return apv In order the subset ind, the regression coefficients, the P-values, the standard P-values.
#' @return res The residuals.
#' @examples 
#' a<-fpval(boston[,14],boston[,1:13],c(1,2,4:6,8:13),13)
fpval<-function(y,x,ind,q=-1,inr=T,xinr=F){
   	n<-length(y)
   	kx<-length(x)/n
   	if(xinr){inr<-F}
   	if(inr){
         		tmpx<-double(n)+1
         		x<-cbind(x,tmpx)
         		ind<-c(ind,kx+1)
         		kx<-kx+1
         		xinr<-TRUE
    	}	
	if(q==-1){q<-kx}
    	ki<-length(ind)
    	apv<-double(2*ki)
    	dim(apv)<-c(ki,2)
    	b<-lm(y~0+x[,ind])
    	res<-as.double(b$res)
    	apv[,2]<-as.double(summary(b)[[4]][1:ki,4])
    	apv[1:ki,1]<-pbeta(apv[1:ki,2],1,q-ki+2)
    	if(xinr){apv[ki,1]<-apv[ki,2]
        		ind[ki]<-0
    	}
   	 beta<-matrix(b$coef,ncol=1)
    	apv<-cbind(ind,beta,apv)
    	list(apv,res)
}
