#' robust stepwise selection of covariates 
#'
#' @param y Dependent variable
#' @param x Covariates
#' @param cn       Parameter fir Huber's psi-function
#' @param alpha  The P-value cut-off
#' @param kmax The maximum number of included covariates
#' @param sig     Scale value of dependent variable
#' @param nu The order statistic of Gaussian covariates used for comparison
#' @param kexk The excluded covariates
#' @param intercept Logical to include intercept
#' @param chkintercept Logical to include or exclude intercept dependent on the P-value
#' @return pv pv[[1]] the included covariates, the P-values; pv[[2]] coefficients of robust linear regression; pv[[3]]  residuals; pv[[4]] scale. 
#' @examples 
#'bostrobreg<-frobstepwise(boston[,14],boston[,1:13],1,0.01,15,intercept=TRUE)
#' bostrobreg[[1]]
#' bostrobreg[[2]]
#' plot(bostrobreg[[3]])
#' bostrobreg[[4]]
frobstepwise<-function(y,x,cn,alpha,nu=1,kmax=0,intercept=TRUE,chkintercept=FALSE,kexk=0,sig=0){
	n<-length(y)
	xx<-x
	if(kmax==0){kmax<-n}
	if(intercept){
		tmpx<-double(n)+1
		dim(tmpx)<-c(n,1)
		xx<-cbind(tmpx,x)
		kmax<-kmax+1
	}
	k<-length(xx[1,])
	kmax1<-kmax+1
	tmpx<-cn*(1:1000)/1000
	cpp<-sum(tmpx^2*dnorm(tmpx))*cn/1000+cn**2*(1-pnorm(cn))
	cpp<-2*cpp
	if(mad(y)==0){stop("MAD=0")}
	if(sig==0){sig<-mad(y)}
	lke<-length(kexk)
	kex<-double(k)
	kex[1:lke]<-kexk
	tmp<-.Fortran(
		"robstepwise",		
		as.double(y),
		as.double(xx),
		as.integer(n),
		as.integer(k),	
		double(n*kmax1),
		double(n),
		double(n),
		double(kmax1),
		integer(k),
		as.double(alpha),
		as.integer(kmax),
		double(2*kmax1),
		double(kmax),
		as.double(cn),
		as.double(cpp),
		as.double(sig),
		double(n),
		double(n),
		double(n),
		as.integer(kmax1),
		as.integer(kex),
		as.logical(intercept),
		as.logical(chkintercept),
		as.double(nu)
	)
	pv<-tmp[[12]]
	dim(pv)<-c(kmax1,2)
	kmax<-tmp[[11]]
	if(kmax==0){
		pv<-c(0,0)
		dim(pv)<-c(1,2)
		res<-NaN
		sig<-NaN
		beta<-NaN	
	}
	else{
		is<-1
		kkmax<-kmax
		if(pv[1,1]==0){is<-2
			kkmax<-kmax-1
		}
		pv[1:kkmax,1]<-pv[is:kmax,1]
		pv[1:kkmax,2]<-pv[is:kmax,2]
		if(kkmax==1){dim(pv)<-c(1,2)
		}
		else{pv<-pv[1:kkmax,]
		}
		sig<-tmp[[16]]
		ind<-pv[,1]
		if(intercept&(!chkintercept)){ind<-ind+1}
		tmpr<-frobreg(y,xx[,ind],cn,cpp=cpp,sig=sig,intercept=intercept)
		res<-tmpr[[2]]
		beta<-tmpr[[1]]
		sig<-tmpr[[3]]
	}
	list(pv,beta,res,sig)
}
