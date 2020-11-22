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
#' @return stpv The in order stepwise P-values, sum of squared residuals and the percentage sum of squared residuals explained.
#' @examples 
#' data(boston)
#' bostint<-fgeninter(boston[,1:13],2)[[1]]
#' a<-f1st(boston[,14],bostint,km=10,sub=T,inr=F,xinr=T)
f1st<-function(y,x,p0=0.01,nu=1,km=0,mx=20,kx=0,sub=F,inr=T,xinr=F){
	if(xinr&(km==1)){stop("only intersect left")}
	n<-length(y)
	x<-matrix(x,nrow=n)
	y<-matrix(y,ncol=1)
	k<-length(x[1,])
	if((!xinr)&inr){tmpx<-double(n)+1
		x<-cbind(x,tmpx)
		x<-matrix(x,nrow=n)
		xinr<-TRUE
	}
	kk<-length(x[1,])
	if(length(kx)==1){
		lkx<-1
		if(kx==0){lkx<-0}
	}
	else{lkx<-length(kx)}
	q<-kk-lkx
	if(xinr){q<-q-1}
	kex<-integer(kk+1)
	if(lkx==1){
		if(kx>0){
			kex[1]<-kx
		}
	}
	if(lkx>1){kex[1:lkx]<-kx}
	if((km>0)&inr){km=km+1}
	p00<-2
	if(km==0){
		p00<-p0
		km<-min(n,kk)
	}
	km1<-km+1
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
		double(k+1)
	)
	kmax<-tmp[[9]]
	if(kmax==0){
		pv<-matrix(c(-2,0,0,0),nrow=1)
		res<-0
		stpv<-0
	}
	else if((kmax==1)&xinr){
		pv<-matrix(c(-1,0,0,0),nrow=1)
		res<-0
		stpv<-0
	}
	else{
		ss01<-1-tmp[[16]][1:kmax]
		stpv<-tmp[[10]]
		stpv<-matrix(stpv,ncol=2)
		stpv<-stpv[1:kmax,]
		stpv<-matrix(stpv,ncol=2)
		minss<-tmp[[15]]
		minss<-minss[1:kmax]
		stpv<-cbind(stpv,minss,ss01)
		stpv<-matrix(stpv,ncol=4)
		ind<-stpv[1:kmax,1]
		if(xinr){ints<-ind[1]
			ind[1:(kmax-1)]<-ind[2:kmax]
			ind[kmax]<-ints
		}
		pv<-fpval(y,x,ind,q,inr=inr,xinr=xinr)
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
				pv<-fpval(y,x,ind,q,inr,xinr)
				li<-length(ind)
				res<-pv[[2]]
				pv<-pv[[1]]
				if(xinr){pv[li,1]<-0}
			}
			else{
				pv<-matrix(c(-1,0,0,0),nrow=1)
				res<-0
				stpv<-0
			}
		}
	}
	if(length(res)==n){
		if(xinr){stpv[1,1]<-0}
	}
	list(pv,res,stpv)
}

