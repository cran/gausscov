#' Calculates all possible subsets and selects those where each included covariate is significant. 
#'
#' It select =TRUE it calls fselect which removes all such subsets which are a subset of some other selected subset. The remaining ones ordered according to the sum of squared residuals
#' 
#'
#' @param y The dependent variable
#' @param x  The covariates
#' @param p0 Cut-off p-value for significance
#' @param select If true use fselect
#' @param intercept If true intercept included
#' @param chkintercept Include intercept depending on p-value
#' @return nv List of subsets with number of covariates and sum of squared residuals
#' @examples 
#' data(redwine)
#' flmmdch(redwine[,12],redwine[,1:11])
#' decode(1874,11)
#' b<-lm(redwine[,12]~redwine[,c(2,5,7,9:11)])
#' summary(b)
#' sum(b$res^2)
#' decode(802,11)
#' b<-lm(redwine[,12]~redwine[,c(2,6,9,10)])
#' summary(b)
#' sum(b$res^2)
flmmdch<-function(y,x,p0=0.01,select=TRUE,intercept=TRUE,chkintercept=FALSE){	
	n<-length(y)
	y<-as.matrix(y)
	x<-as.matrix(x)
	k<-length(x[1,])
	nv<-integer(2^(k+1))
	dim(nv)=c(2^k,2)
	if(intercept){
		if(!chkintercept){y<-y-mean(y)
			for(j in 1:k){
				x[,j]<-x[,j]-mean(x[,j])
			}
		}
		else{tmpx<-double(n)+1
			x<-cbind(tmpx,x)
			k<-k+1
		}
	}
	nv<-integer(2^(k+1))
	dim(nv)=c(2^k,2)
	tmp<-.Fortran(
		"lmmdch",
		as.double(y),
		as.double(x),
		as.integer(n),
		as.integer(k),	
		double(n*k),
		double(n*k),
		double(n),
		double(n),
		double(k),
		double(k),
		double(k),
		double(k^2),
		integer(k),
		as.logical(intercept),
		double(2^k),
		as.integer(nv),
		double(2^k),
		as.double(p0)
	)
	ss<-tmp[[17]]
	nv<-tmp[[16]]
	dim(nv)<-c(2^k,2)
	inv<-(1:2^k)[nv[,2]>0]
	nv<-nv[inv,]
	ss<-ss[inv]
	inv<-1:length(inv)
	ind<-rank(ss)
	inv[ind]<-inv
	nv<-nv[inv,]
	ss<-ss[inv]
#
#	select models
#
	nv<-cbind(nv,ss)
	if(select){
		ind<-fselect(nv,k)[[1]]
		nv<-nv[ind,]
	}
	list(nv)
}

