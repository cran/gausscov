#' Selection of covariates with given excluded covariates, called in f3st.R 
#'
#' @param y Dependent variable.
#' @param x Covariates.
#' @param covch Sum of squared residuals and selected covariates
#' @param ind The excluded covariates
#' @param m Number of iterations
#' @param kexmx The maximum number of covariates in an approximation.
#' @param p0 The cut-off P-value.
#' @param nu The order statistic of Gaussian covariates used for comparison.
#' @param kmn The minimum number of included covariates irrespective of cut-off P-value.
#' @param kmx The maximum number of included covariates irrespective of cut-off P-value.
#' @param mx  The maximum number covariates for an all subset search.
#' @param lm  The maximum number of approximations.
#' @param kex The excluded covariates.
#' @param sub Logical, if TRUE best subset selected.
#' @param inr Logical, if TRUE include intercept if not present.
#' @param xinr Logical, if TRUE intercept already included 
#' @qq   The number of covariates to choose from. If qq=0 the number of covariates of x is used.
#' @lm0 The current number of approximations.
#' @return ind1 The excluded covariates
#' @return covch The sum of squared residuals and the selected covariates ordered in increasing size of sum of squared residuals.
#' @returm lm0 The current number of approximations.
#' @examples 
#' data(leukemia)
#' covch<-c(2.023725,1182,1219,2888,0)
#' covch<-matrix(covch,nrow=1,ncol=5)
#' ind<-c(1182,1219,2888)
#' ind<-matrix(ind,nrow=3,ncol=1)
#' m<-1
#' a<-f3sti(leukemia[[1]],leukemia[[2]],covch,ind,m)
f3sti<-function(y,x,covch,ind,m,kexmx=100,p0=0.01,nu=1,kmn=0,kmx=0,mx=21,lm=1000,kex=0,sub=T,inr=T,xinr=F,qq=0,lm0=0){
	ind<-matrix(ind,ncol=m)
	ni<-length(ind[,1])
	ind1<-integer(m+1)
	ind1<-matrix(ind1,ncol=m+1,nrow=1)
	kexx<-integer(m)
	kexx<-matrix(kexx,ncol=m)
	i<-1
	while(i<=ni){
		kexi<-sort((1:m)[ind[i,]>0])
		lkx<-length(kexx[,1])
		kex<-ind[i,kexi]
		rp<-1
		iu<-1
		while(iu<=lkx){
			tmp<-max(abs(kexx[iu,]-kex))
			if(tmp==0){rp<-0
				iu<-lkx
			}
			iu<-iu+1
		}
		if((rp>0)&lm0<lm){
			lkex<-length(kex)
			covch0<-double(kexmx)
			tt<-system.time(a<-f1st(y,x,p0=p0,kmn=kmn,sub=sub,kex=kex,inr=inr,xinr=xinr))[[1]]
			if(a[[1]][1,1]>0){
				la<-length(a[[1]][,1])-1
				cont<-TRUE
				if(la>=(kexmx-1)){
					cont<-FALSE 
					warning("Number of covariates exceeds limit kexmx")
				}
				if(cont&(la>0)){
					ssa2<-sum(a[[2]]^2)
					ssmin<-min(abs(covch[,1]-ssa2))
					if(ssmin>0){
						covch0[1]<-ssa2
						covch0[2:(la+1)]<-a[[1]][1:la,1]
						kexx<-rbind(kexx,kex)
						covch<-rbind(covch,covch0)
						lm0<-lm0+1
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
		}
		i<-i+1
	}
	if(sum(ind1)>0.5){
		il<-length(ind1[,1])
		ind1<-ind1[2:il,]
		covch<-matrix(covch,ncol=kexmx)
	}
	li<-length(covch[,1])
	covch<-matrix(covch,nrow=li,ncol=kexmx)
	ind0<-rank(covch[,1],ties.method="first")
	lli<-1:li
	lli[ind0]<-lli
	covch<-covch[lli,]
	list(ind1,covch,lm0)
}
