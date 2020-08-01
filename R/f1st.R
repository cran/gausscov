#' Stepwise selection of covariates 
#'
#' @param y Dependent variable
#' @param x Covariates
#' @param p0 The P-value cut-off
#' @param nu The order statistic of Gaussian covariates used for comparison
#' @param km The maximum number of included covariates
#' @param mx  The maximum number covariates for an all subset search.
#' @param kx The excluded covariates
#' @param sub Logical, if T best subset selected.
#' @param inr Logical, if T include intercept if not present
#' @param xinr Logical, if T intercept already present 
#' @return pv The in order selected covariates, the regression coefficients, the P-values, the standard P-values.
#' @return res The residuals.
#' @return stpv The stepwise P-values and the sum of squared residuals
#' @examples 
#' data(boston)
#' bostint<-fgeninter(boston[,1:13],2)[[1]]
#' a<-f1st(boston[,14],bostint,km=10,sub=TRUE,inr=FALSE,xinr=TRUE)
f1st<-function(y,x,p0=0.01,nu=1,km=0,mx=20,kx=0,sub=F,inr=T,xinr=F){
	if(xinr&(km==1)){stop("only intersect left")}
	n<-length(y)
	x<-matrix(x,nrow=n)
	y<-matrix(y,ncol=1)
	k<-length(x[1,])
	q<-k
	if(xinr){q<-q-1}
	if((!xinr)&inr){tmpx<-double(n)+1
		x<-cbind(x,tmpx)
		x<-matrix(x,nrow=n)
		xinr<-TRUE
	}
	kk<-length(x[1,])
	kex<-integer(kk+1)
	lke<-length(kx)
	if(lke==1){
		if(kx>0){
			kex[1]<-kx
		}
	}
	if(lke>1){kex[1:lke]<-kx}
	if((km>0)&inr){km=km+1}
	p00<-2
	if(km==0){
		p00<-p0
		km<-min(n,kk)
	}
	km1<-km+1
	pv<-double(3*kk)
	dim(pv)<-c(kk,3)
	tmp<-.Fortran(
		"fstepwise",
		as.double(y),
		as.double(x),
		as.integer(n),
		as.integer(kk),	
		double(n),
		double(n),
		integer(kk+1),
		as.double(p00),
		as.integer(km),
		double(2*km1),
		as.integer(km1),
		as.integer(kex),
		as.logical(xinr),
		as.double(nu),
		double(km1),
		double(k)
	)
	kmax<-tmp[[9]]
	if(kmax==0){
		pv<-matrix(c(-2,0,0,0),nrow=1)
		res<-0
		sig<-0
		stpv<-0
	}
	else if((kmax==1)&xinr){
		pv<-matrix(c(-1,0,0,0),nrow=1)
		res<-0
		sig<-0
		stpv<-0
	}
	else{
		ss01<-tmp[[16]][1:kmax]
		stpv<-tmp[[10]]
		stpv<-matrix(stpv,ncol=2)
		stpv<-stpv[1:kmax,]
		minss<-tmp[[15]]
		minss<-minss[1:kmax]
		stpv<-cbind(stpv,minss,ss01)
		stpv<-matrix(stpv,ncol=4)
		ind<-stpv[1:kmax,1]
		ind<-sort(ind)
		pv<-fpval(y,x,ind,q,xinr)
		res<-pv[[2]]
		pv<-pv[[1]]
		li<-length(ind)
		if(xinr){pv[li,1]<-0}
		if(sub&(kmax>=2)){
			if(kmax>mx){stop("kmax too large (> mx) for all subset search")}
			sbsts<-fmch(y,x,p0=p0,q=q,ind=ind,inr=inr,xinr=xinr)[[1]]
			if(sbsts[1,1]>1){
				tmv<-decode(sbsts[1,1],k)[[1]]
				ind<-ind[tmv]
				if(xinr){ind<-c(ind,kk)}
				pv<-fpval(y,x,ind,q,xinr=xinr)
				li<-length(ind)
				res<-pv[[2]]
				pv<-pv[[1]]
				if(xinr){pv[li,1]<-0}
			}
		}
	}
	list(pv,res,stpv)
}

