#' Calculates the regression coefficients, the P-values and the standard P-values for the chosen subset ind
#'
#' @param y The dependent variable
#' @param x  The covariates
#' @param ind The indices of the subset whose P-values are required
#' @param q   The total  number of covariates from which ind was selected
#' @param xinr Logical If TRUE intercept included
#' @return apv In order the subset ind, the regression coefficients, the P-values, the standard P-values.
#' @return res The residuals.
#' @examples 
#' a<-fpval(boston[,14],boston[,1:13],c(1,2,4:6,8:13),13)
fpval<-function(y,x,ind,q,xinr=F){
	k<-length(ind)
	apv<-double(2*k)
	dim(apv)<-c(k,2)
	b<-lm(y~0+x[,ind])
	res<-as.double(b$res)
	apv[,2]<-as.double(summary(b)[[4]][1:k,4])
	apv[1:k,1]<-pbeta(apv[1:k,2],1,q-k+2)
	if(xinr){apv[k,1]<-apv[k,2]}
	beta<-matrix(b$coef,ncol=1)
	apv<-cbind(ind,beta,apv)
	list(apv,res)		
}
