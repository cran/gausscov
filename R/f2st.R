#' Repeated stepwise selection of covariates 
#'
#' @param y Dependent variable.
#' @param x Covariates.
#' @param p0  The cut-off P-value.
#' @param kmn The minimum number of included covariates irrespective of cut-off P-value.
#' @param kmx The maximum number of included covariates irrespective of cut-off P-value.
#' @param kex The excluded covariates.
#' @param mx  The maximum number covariates for an all subset search.
#' @param lm  The maximum number of linear approximations. 
#' @param sub Logical, if T select best subset.
#' @param inr Logical if T include an intercept .
#' @param xinr Logical if T intercept already included.
#' @param qq   The number of covariates to choose from. If qq=0 the number of covariates of x is used.
#' @return pv   In order the linear approximation,  the included covariates, the regression coefficient values, the Gaussian P-values, the standard P-values  and the proportional reduction in the sum of squared residuals due to this covariate.
#' @examples 
#' data(boston)
#' bostint<-fgeninter(boston[,1:13],2)[[1]]
#' a<-f2st(boston[,14],bostint,lm=3,sub=FALSE)
f2st<-function(y,x,p0=0.01,kmn=0,kmx=0,kex=0,mx=21,lm=9^9,sub=T,inr=T,xinr=F,qq=0){
	kxx<-kex
	kmmx<-kmx
	kmnn<-kmn
	n<-length(y)
	if((!xinr)&inr){
		tmpx<-double(n)+1
		x<-cbind(x,tmpx)
		xinr<-TRUE
	}
	q<-length(x[1,])
	if(xinr){q<-q-1}
	kex<-0
	kv<-1
	mn<-0
	qqq<-qq
	pvv<- matrix(c(0,0,0,0,0),nrow=1)
	while(kv>0.5){
		tmp<-f1st(y,x,p0=p0,kmn=kmnn,kmx=kmmx,kex=kxx,mx=mx,sub=sub,inr=inr,xinr=xinr,qq=qq)[[1]]
		if(tmp[1,1]<=0){kv<-0}
		if(kv>0.5){
			lv<-length(tmp[,1])
			llv<-(1:lv)[tmp[,1]>0]
			mn<-mn+1
			mnt<-integer(lv)+mn
			mnt<-matrix(mnt,ncol=1)
			if(mn==1){pv<-tmp
				pvv<-cbind(mnt,pv)
				if(length(kxx)==1){
					if(kxx==0){
						kxx<-tmp[llv,1]
					}
					else{
						kxx<-c(kxx,tmp[llv,1])
					}
				}
				else{kxx<-c(kxx,tmp[llv,1])}
				mnn<-mnt
			}
			else{mnn<-rbind(mnn,mnt)
				pv<-rbind(pv,tmp)
				pvv<-cbind(mnn,pv)
				kxx<-c(kxx,tmp[llv,1])
			}
			kv<-1
			if(mn>=lm){kv<-0}
			if(length(kxx)>=q){kv<-0}
		}
	}
	list(pvv)
}
