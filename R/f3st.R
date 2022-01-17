#' Stepwise exclusion and selection of covariates 
#'
#' @param y Dependent variable.
#' @param x Covariates.
#' @param m The number of iterations
#' @param kexmx The maximum number of excluded covariates
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
#' @return covch The sum of squared residuals and the selected covariates ordered in increasing size of sum of squared residuals.
#' @return lai The length of covch
#' @examples 
#' data(leukemia_y)
#' data(leukemia_x)
#' a<-f3st(ly.original,lx.original,m=2,kexmx=5,kmn=5,sub=TRUE)
f3st<-function(y,x,m,kexmx=100,p0=0.01,nu=1,kmn=0,kmx=0,mx=21,kex=0,sub=T,inr=T,xinr=F,qq=0){	
	covch<-double(kexmx)
	covch<-matrix(covch,nrow=1,ncol=kexmx)
	a<-f1st(y,x,p0=p0,kmn=kmn,sub=T)
	if(a[[1]][1,1]>0){
		ali<-length(a[[1]][,1])-1
		if(ali>0){
			covch[1]<-sum(a[[2]]^2)
			a<-a[[1]]
			covch[2:(ali+1)]<-a[1:ali,1]
			ind<-a[1:ali,1]
			ind<-matrix(ind,ncol=1)
			li<-length(ind[,1])
			for(j in 1:m){
				a<-f3sti(y,x,covch,ind,j,kexmx=kexmx,p0=p0,nu=nu,kmn=kmn,kmx=kmx,mx=mx,kex=kex,sub=sub,inr=inr,xinr=xinr,qq=qq)
				if(sum(a[[1]])>0){ind<-a[[1]]}
				covch<-a[[2]]
			}
			li<-length(covch)/kexmx
			covch<-matrix(covch,nrow=li,ncol=kexmx)
			if(li>1){
				ind<-rank(covch[,1],ties.method="first")
				lli<-1:li
				lli[ind]<-lli
				covch<-covch[lli,]
				ind<-diff(covch[,1])
				li<-length(ind)
				aind<-(1:li)[ind>1e-06]
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
	list(covch,lai)
}


