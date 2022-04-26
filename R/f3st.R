#' Stepwise exclusion and selection of covariates 
#'a<-cv.glmnet(x,yy))a<-cv.glmnet(x,yy))g
#' @param y Dependent variable.
#' @param x Covariates.
#' @param m The number of iterations
#' @param kexmx The maximum number of covariates in an approximation
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
#' @return covch The sum of squared residuals and the selected covariates ordered in increasing size of sum of squared residuals.
#' @return lai The number of rows of covch
#' @examples 
#' data(leukemia)
#' a<-f3st(leukemia[[1]],leukemia[[2]],m=1,kexmx=5,kmn=5)
f3st<-function(y,x,m,kexmx=100,p0=0.01,nu=1,kmn=0,kmx=0,mx=21,lm=1000,kex=0,sub=T,inr=T,xinr=F,qq=0){
	kex0<-0
	if(sum(kex)>0){kex0<-sort(kex)
		x<-x[,-kex0]
		kex<-0
	}
	covch<-double(kexmx)
	covch<-matrix(covch,nrow=1,ncol=kexmx)
	a<-f1st(y,x,p0=p0,nu=nu,kmn=kmn,kmx=kmx,mx=mx,kex=kex,sub=sub,inr=inr,xinr=xinr,qq=qq)
	if(a[[1]][1,1]>0){
		ali<-length(a[[1]][,1])-1
		if(ali>0){
			covch[1]<-sum(a[[2]]^2)
			a<-a[[1]]
			covch[2:(ali+1)]<-a[1:ali,1]
			ind<-a[1:ali,1]
			ind<-sort(ind)
			ind<-matrix(ind,ncol=1)
			li<-length(ind[,1])
			lm0<-1
			j<-1
			while(j<=m){
				a<-f3sti(y,x,covch,ind,j,kexmx=kexmx,p0=p0,nu=nu,kmn=kmn,kmx=kmx,mx=mx,lm=lm,kex=kex,sub=sub,inr=inr,xinr=xinr,qq=qq,lm0=lm0)
				if(sum(a[[1]])>0){ind<-a[[1]]}
				covch<-a[[2]]
				lm0<-a[[3]]
				if(lm0==lm){j<-m}
				j<-j+1
				if(sum(a[[1]])==0){j<-m+1}
			}
			li<-length(covch)/kexmx
			covch<-matrix(covch,nrow=li,ncol=kexmx)
			lc<-length(covch)/kexmx
			if(li>1){
				ind<-rank(covch[,1],ties.method="first")
				lli<-1:li
				lli[ind]<-lli
				covch<-covch[lli,]
				ind<-diff(covch[,1])
				li<-length(ind)
				aind<-(1:li)[ind>1e-10]+1
				aind<-c(1,aind)
				covch<-covch[aind,]
				lc<-length(covch)/kexmx
				covch<-matrix(covch,nrow=lc,ncol=kexmx)
				laind<-length(aind)
				j<-1
				while(j<=kexmx){
					if(lc==1){
						jj<-covch[1,j]
						if(jj<0.5){jm<-j-1
							j<-kexmx+1
						}
					}
					else{
						jj<-max(covch[1:lc,j])
						if(jj<0.5){jm<-j-1
							j<-kexmx+1
						}
					}
					j<-j+1
				}
				covch<-covch[,1:jm]
				covch<-matrix(covch,nrow=lc,ncol=jm)
				lai<-jm
				
			}
			else{
				li<-(1:kexmx)[covch[1,1:kexmx]>0]
				covch<-covch[1,li]
				lai<-length(li)
			}

		}
	}
	else{covch<-c(-1,0)
		lai<-0
	}
	if((lai>0)&(sum(kex0)>0)){
		lc<-length(covch[,1])
		for(j in kex0){
			for(i in 1:lc){
				for(ij in 2:lai){
					if(covch[i,ij]>=j){covch[i,ij]<-covch[i,ij]+1}
				}
			}
		}
	}
	list(covch,lai)
}


