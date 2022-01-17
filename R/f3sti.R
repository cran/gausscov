#' Selection of covariates with given excluded covariates, called in f3st.R 
#'
#' @param y Dependent variable.
#' @param x Covariates.
#' @param covch Sum of squared residuals and selected covariates
#' @param ind The excluded covariates
#' @param m Number of iterations
#' @param kexmx mmaximum number of excluded covariates.
#' @param p0 The cut-off P-value.
#' @param nu The order statistic of Gaussian covariates used for comparison.
#' @param kmn The minimum number of included covariates irrespective of cut-off P-value.
#' @param kmx The maximum number of included covariates irrespective of cut-off P-value.
#' @param mx  The maximum number covariates for an all subset search.
#' @param kex The excluded covariates.
#' @param sub Logical, if TRUE best subset selected.
#' @param inr Logical, if TRUE include intercept if not present.
#' @param xinr Logical, if TRUE intercept already included 
#' @qq   The number of covariates to choose from. If qq=0 the number of covariates of x is used.
#' @return ind1 The excluded covariates
#' @return covch The sum of squared residuals and the selected covariates ordered in increasing size of sum of squared residuals.
#' @examples 
#' data(leukemia_y)
#' data(leukemia_x)
#' covch=c(2.023725,1182,1219,2888,0)
#' ind<-c(1182,1219,2888,0,0,0,0,0,0)
#' ind<-matrix(ind,ncol=3)
#' m<-3
#' a<-f3sti(ly.original,lx.original,covch,ind,m)
f3sti<-function(y,x,covch,ind,m,kexmx=100,p0=0.01,nu=1,kmn=0,kmx=0,mx=21,kex=0,sub=T,inr=T,xinr=F,qq=0){
	ind<-matrix(ind,ncol=m)
	ni<-length(ind[,1])
	ind1<-integer(m+1)
	ind1<-matrix(ind1,ncol=m+1,nrow=1)
	for(i in 1:ni){
		kexi<-(1:m)[ind[i,]>0]
		kex<-ind[i,kexi]
		lkex<-length(kex)
		covch0<-double(kexmx)
		a<-f1st(y,x,p0=p0,kmn=kmn,sub=sub,kex=kex)
		if(a[[1]][1,1]>0){
			la<-length(a[[1]][,1])-1
			if(la>0){
				covch0[1]<-sum(a[[2]]^2)
				covch0[2:(la+1)]<-a[[1]][1:la,1]
				covch<-rbind(covch,covch0)
				for(iu in 1:la){
					ind0<-integer(m+1)
					ind0<-matrix(ind0,nrow=1,ncol=m+1)
					ind0[1:lkex]<-kex
					ind0[(1+lkex)]<-a[[1]][iu,1]
					ind1<-rbind(ind1,ind0)
				}
			}
		}
	}
	if(sum(ind1)>0.5){
		il<-length(ind1[,1])
		ind1<-ind1[2:il,]
		covch<-matrix(covch,ncol=kexmx)
	}
	list(ind1,covch)
}
