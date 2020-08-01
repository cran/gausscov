#' Robust stepwise selection of covariates 
#'
#' @param y         Dependent variable
#' @param x         Covariates
#' @param cn        Parameter for Huber's psi-function
#' @param p0       The P-value cut-off
#' @param sg       Scale value of residuals
#' @param nu       The order for calculating the P-value
#' @param km      The maximum number of included covariates
#' @param mx      The maximum number of included covariates if the option subset =TRUE is used.
#' @param kx       The excluded covariates
#' @param sb     Logical If TRUE best subset selected 
#' @param of Logical if TRUE to include intercept 
#' @param xof Logical if TRUE intercept already included 
#' @return pv    In order the subset ind, the regression coefficients, the P-values, the standard P-values.
#' @return res   The residuals
#' @return stpv  The stepwise regression results: covariate, P-value and scale. 
#' @examples 
#' data(boston)
#' a<-frst(boston[,14],boston[,1:13],km=10,sb=T)
frst<-function(y,x,cn=1,p0=0.01,sg=0,nu=1,km=0,mx=20,kx=0,sb=F,of=T,xof=F){
	if(mad(y)==0){stop("MAD=0")}
	if(sg==0){sg<-mad(y)}
	tmpx<-cn*(1:1000)/1000
	cpp<-sum(tmpx^2*dnorm(tmpx))*cn/1000+cn**2*(1-pnorm(cn))
	cpp<-2*cpp
	n<-length(y)
	x<-matrix(x,nrow=n)
	k<-length(x)/n
	q<-k
	if(xof){q<-q-1}
	if((!xof)&of){tmpx<-double(n)+1
		x<-cbind(x,tmpx)
		x<-matrix(x,nrow=n)
		xof<-T
	}
	k<-length(x[1,])
	kex<-integer(k+1)
	lke<-length(kx)
	if(lke==1){
		if(kx>0){
			kex[1]<-kx
		}
	}
	if(lke>1){kex[1:lke]<-kx}
	if((km>0)&of){km=km+1}
	p00<-2
	if(km==0){
		p00<-p0
		km<-min(n,k)
	}
	km1<-km+1
	pv<-double(3*k)
	dim(pv)<-c(k,3)
	tmp<-.Fortran(
		"robstepwise",		
		as.double(y),
		as.double(x),
		as.integer(n),
		as.integer(k),	
		double(n*km1),
		double(n),
		double(n),
		double(km1),
		integer(k),
		as.double(p00),
		as.integer(km),
		double(2*km1),
		double(km1),
		as.double(cn),
		as.double(cpp),
		as.double(sg),
		double(n),
		double(n),
		double(n),
		as.integer(km1),
		as.integer(kex),
		as.logical(xof),
		as.double(nu),
		double(km1)
	)
	kmax<-tmp[[11]]
	if(kmax==0){
		pv<-matrix(c(-2,0,0,0),nrow=1)
		res<-0
		sig<-0
		stpv<-0
	}
	else if((kmax==1)&xof){
		pv<-matrix(c(-1,0,0,0),nrow=1)
		res<-0
		sig<-0
		stpv<-0
	}
	else{
		ss01<-tmp[[24]][1:kmax]
		stpv<-tmp[[12]]
		stpv<-matrix(stpv,ncol=2)
		stpv<-stpv[1:kmax,]
		sg<-tmp[[16]]
		stpv<-cbind(stpv,ss01)
		stpv<-matrix(stpv,ncol=3)
		ind<-stpv[1:kmax,1]
		ind<-sort(ind)
		pv<-fpval(y,x,ind,q,xof)
		res<-pv[[2]]
		pv<-pv[[1]]
		li<-length(ind)
		if(xof){pv[li,1]<-0}
		if(sb&(kmax>=2)){
			if(kmax>mx){stop("kmax too large (> mx) for all subset search")}
			sbsts<-fmch(y,x,p0=p0,q=q,ind=ind,inr=of,xinr=xof)[[1]]
			if(sbsts[1,1]>1){
				tmv<-decode(sbsts[1,1],k)[[1]]
				ind<-ind[tmv]
				if(xof){ind<-c(ind,k)}
				pv<-fpval(y,x,ind,q,xinr=xof)
				li<-length(ind)
				res<-pv[[2]]
				pv<-pv[[1]]
				if(xof){pv[li,1]<-0}
			}
		}
	}
	list(pv,res,stpv)
}

