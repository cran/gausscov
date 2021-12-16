#' Stepwise selection of covariates 
#'
#' @param y Dependent variable.
#' @param x Covariates.
#' @param p0 The cut-off P-value.
#' @param nu The order statistic of Gaussian covariates used for comparison.
#' @param kmn The minimum number of included covariates irrespective of cut-off P-value.
#' @param kmx The maximum number of included covariates irrespective of cut-off P-value.
#' @param mx  The maximum number covariates for an all subset search.
#' @param kex The excluded covariates.
#' @param sub Logical, if TRUE best subset selected.
#' @param inr Logical, if TRUE include intercept if not present.
#' @param xinr Logical, if TRUE intercept already included 
#' @param qq   The number of covariates to choose from. If qq=0 the number of covariates of x is used.
#' @return pv In order the included covariates, the regression coefficient values, the Gaussian P-values, the standard P-values  and the proportional reduction in the sum of squared residuals due to this covariate
#' @return res The residuals.
#' @return stpv The in order stepwise P-values, sum of squared residuals and the proportional  reduction in the sum of squared residuals due to this covariate.
#' @examples 
#' data(boston)
#' bostint<-fgeninter(boston[,1:13],2)[[1]]
#' a<-f1st(boston[,14],bostint,kmn=10,sub=TRUE)
f1st<-function(y,x,p0=0.01,nu=1,kmn=0,kmx=0,mx=21,kex=0,sub=T,inr=T,xinr=F,qq=0){
	if(xinr&(kmn==1)){stop("only intersect left")}
	n<-length(y)
	x<-matrix(x,nrow=n)
	y<-matrix(y,ncol=1)
	kex<-matrix(kex,nrow=1)
	lkx<-length(kex)
	if(!xinr){
		if(inr){
			tmpx<-double(n)+1
			x<-cbind(x,tmpx)
			x<-matrix(x,nrow=n)
			xinr<-TRUE
			inr<-FALSE
		}
	}
	k<-length(x[1,])
	if(xinr){
		if(kex[1]==0){q<-k-1}
		else{q<-k-1-lkx}
	}
	else{
		if(kex[1]==0){q<-k}
		else{q<-k-lkx}
	}
	if(kex[1]==0){
		kex<-integer(k)
	}
	else{tmp<-kex
		kex<-integer(k)
		kex[1:lkx]<-tmp
	}
	if((kmn>0)&xinr){kmn=kmn+1}
	if(kmx==0){kmx<-min(n,k)}
	tmp<-.Fortran(
		"fstepwise",
		as.double(y),
		as.double(x),
		as.integer(n),
		as.integer(k),	
		double(n),
		double(n),
		integer(k),
		as.double(p0),
		as.integer(kmx),
		double(2*(k+1)),
		as.integer(kex),
		as.logical(xinr),
		as.double(nu),
		double(k),
		double(k),
		as.integer(qq),
		as.integer(kmn)
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
		ss01<-1-tmp[[15]][1:kmax]
		stpv<-tmp[[10]]
		stpv<-matrix(stpv,ncol=2)
		stpv<-stpv[1:kmax,]
		stpv<-matrix(stpv,ncol=2)
		minss<-tmp[[14]]
		minss<-minss[1:kmax]
		stpv<-cbind(stpv,minss,ss01)
		stpv<-matrix(stpv,ncol=4)
		ind<-stpv[1:kmax,1]
		if(xinr){ints<-ind[1]
			ind[1:(kmax-1)]<-ind[2:kmax]
			ind[kmax]<-ints
		}
		if(qq==0){
			pv<-fpval(y,x,ind,q=q,nu=nu,inr=inr,xinr=xinr)
		}
		else{
			pv<-fpval(y,x,ind,q=qq,nu=nu,inr=inr,xinr=xinr)
		}
		res<-pv[[2]]
		pv<-pv[[1]]
		li<-length(ind)
		if(xinr){pv[li,1]<-0}
		ssb<--1
		if((kmax>mx)&(sub==TRUE)){sub<-FALSE
			warning("kmax too large (> mx) for all subset search")
		}
		if(sub&(kmax>=2)){
			if(qq==0){
				sbsts<-fasb(y,x,p0=p0,q=q,ind=ind,inr=inr,xinr=xinr,)[[1]]
			}
			else{
				sbsts<-fasb(y,x,p0=p0,q=qq,ind=ind,inr=inr,xinr=xinr,)[[1]]
			}
			sbsts<-matrix(sbsts,nrow=3)
			kmm<-length(ind)
			if(xinr){kmm<-kmm-1}
			if(sbsts[1,1]>0){
				nss<-sbsts[1,1]
				tmv<-decode(nss,kmm)[[1]]
				ind<-ind[tmv]
				if(xinr){if(max(ind)<k){ind<-c(ind,k)}}
				if(qq==0){
					a<-fpval(y,x,ind,q=q,nu=nu,inr=inr,xinr=xinr)
				}
				else{
					a<-fpval(y,x,ind,q=qq,nu=nu,inr=inr,xinr=xinr)
				}
				li<-length(ind)
				res<-a[[2]]
				pv<-a[[1]]
				if(xinr){pv[li,1]<-0}
			}
			else{
				pv<-matrix(c(-1,0,0,0,0),nrow=1)
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

