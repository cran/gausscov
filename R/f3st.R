#' Stepwise exclusion and selection of covariates 
#'
#' @param y Dependent variable.
#' @param x Covariates.
#' @param m The number of iterations
#' @param kmn The minimum number of included covariates irrespective of cut-off P-value.
#' @param p0 The cut-off P-value.
#' @param kmx The maximum number of included covariates irrespective of cut-off P-value.
#' @param mx  The maximum number covariates for an all subset search.
#' @param lm  The maximum number of approximations.
#' @param kex The excluded covariates.
#' @param sub Logical, if TRUE best subset selected.
#' @param inr Logical, if TRUE include intercept if not present.
#' @param xinr Logical, if TRUE intercept already included 
#' @qq   The number of covariates to choose from. If qq=0 the number of covariates of x is used.
#' @param kexmx The maximum number of covariates in an approximation
#' @return covch The sum of squared residuals and the selected covariates ordered in increasing size of sum of squared residuals.
#' @return lai The number of rows of covch
#' @examples 
#' data(leukemia)
#' a<-f3st(leukemia[[1]],leukemia[[2]],m=1,kmn=5,kexmx=5)
f3st<-function(y,x,m,kmn=10,p0=0.01,kmx=0,mx=21,lm=100,kex=0,sub=T,inr=T,xinr=F,qq=0,kexmx=100){
	kex0<-0
	if(sum(kex)>0){kex0<-sort(kex)
		x<-x[,-kex0]
		kex<-0
	}
	covch<-double(kexmx)
	covch<-matrix(covch,nrow=1,ncol=kexmx)
	a<-f1st(y,x,p0=p0,kmn=kmn,kmx=kmx,mx=mx,kex=kex,sub=sub,inr=inr,xinr=xinr,qq=qq)
	if(a[[1]][1,1]>0){
		ali<-length(a[[1]][,1])
		ind<-(1:ali)[a[[1]][,1]>0]
		ind<-a[[1]][ind,1]
		ali<-length(ind)
		covch[1,1]<-sum(a[[2]]^2)
		a<-a[[1]]
		ind<-sort(ind)
		covch[1,2:(ali+1)]<-ind
		li<-length(ind)
		lm0<-1
		lai<-1
		if(m>0){
			a<-f3sti(y,x,covch,ind,m,kexmx=kexmx,p0=p0,kmn=kmn,kmx=kmx,mx=mx,lm=lm,kex=kex,sub=sub,inr=inr,xinr=xinr,qq=qq,lm0=lm0)
		covch<-a[[1]]
		lai<-length(covch[,1])
		}
	}
	else{lai<-0
		covch<-0
	}
	tmp<-apply(covch,2,max)
	len<-length(tmp)
	ind<-(1:len)[tmp>0]
	covch<-covch[,ind]
	list(covch,lai)
}


