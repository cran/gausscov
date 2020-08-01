#' Repeated stepwise selection of covariates 
#'
#' @param y Dependent variable
#' @param x Covariates
#' @param p0  The P-value cut-off
#' @param nu The order statistic of Gaussian covariates used for comparison
#' @param km The maximum specified number of included covariates at each stage
#' @param lm  The maximum number of linear approximations 
#' @param mx  The maximum number covariates for an all subset search.
#' @param kx The excluded covariates
#' @param sub Logical, if TRUE best subset selected.
#' @param inr Logical if TRUE to include an intercept 
#' @param xinr Logical if TRUE intercept already included
#' @return pv   In order, the linear approximation, the included covariates, the coefficient values, the P-values and the standard P-values
#' @examples 
#' data(boston)
#' bostint<-fgeninter(boston[,1:13],2)[[1]]
#' a<-f2st(boston[,14],bostint,lm=3,sub=T)
f2st<-function(y,x,p0=0.01,nu=1,km=0,mx=20,kx=0,lm=9^9,sub=F,inr=T,xinr=F){
	kk<-length(x[1,])+1
	n<-length(y)
	q<-length(x[1,])
	kx<-0
	kv<-1
	mn<-0
	while(kv>0.5){
		tmp<-f1st(y,x,p0,nu,km,mx,kx,sub,inr,xinr)[[1]]
		if((!xinr)&inr){xinr<-T}
		if(tmp[1,1]<=0){kv<-0}
		if(kv>0.5){
			lv<-length(tmp[,1])
			llv<-(1:lv)[tmp[,1]>0]
			mn<-mn+1
			mnt<-integer(lv)+mn
			mnt<-matrix(mnt,ncol=1)
			if(mn==1){pv<-tmp
				pvv<-cbind(mnt,pv)
				kx<-tmp[llv,1]
				mnn<-mnt
			}
			else{mnn<-rbind(mnn,mnt)
				pv<-rbind(pv,tmp)
				pvv<-cbind(mnn,pv)
				kx<-c(kx,tmp[llv,1])
			}
			kv<-1
			if(mn>=lm){kv<-0}
			if(length(kx)>=q){kv<-0}
		}
	}
	list(pvv)
}
