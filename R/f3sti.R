' Selection of covariates with given excluded covariates, called in f3st.R 
#'
#' @param y Dependent variable.
#' @param x Covariates.
#' @param covch Sum of squared residuals and selected covariates
#' @param ind The excluded covariates
#' @param m Number of iterations
#' @param p0 The cut-off P-value.
#' @param kmn The minimum number of included covariates irrespective of cut-off P-value.
#' @param kmx The maximum number of included covariates irrespective of cut-off P-value.
#' @param kex The excluded covariates.
#' @param mx  The maximum number covariates for an all subset search.
#' @param sub Logical, if TRUE best subset selected.
#' @param inr Logical, if TRUE include intercept if not present.
#' @param xinr Logical, if TRUE intercept already included 
#' @qq   The number of covariates to choose from. If qq=-1 the number of covariates of x is used.
#' @param kexmx The maximum number of covariates in an approximation
#' @return covch The sum of squared residuals and the selected covariates ordered in increasing size of sum of squared residuals.
#' @examples 
#' data(leukemia)
#' covch<-c(2.023725,1182,1219,2888,0)
#' covch<-matrix(covch,nrow=1,ncol=5)
#' ind<-c(1182,1219,2888)
#' ind<-matrix(ind,nrow=3,ncol=1)
#' m<-1
#' a<-f3sti(leukemia[[1]],leukemia[[2]],covch,ind,m,kexmx=5)
f3sti<-function(y,x,covch,ind,m,p0=0.01,kmn=0,kmx=0,kex=0,mx=21,sub=T,inr=F,xinr=F,qq=-1,kexmx=100){
	kexx<-integer(m)
	kexx<-matrix(ind,ncol=1)
	if(max(kexx)>0){
		lkex<-length(kex)
		kexx[1:lkex]<-kex
	}
	mm<-1
	while(mm <= m){
			ind<-kexx
			ind<-matrix(ind,ncol=mm)
			ni<-length(ind[,1])
			kex0<-integer(mm+1)
			for(i in 1:ni){
				kex<-ind[i,]
				a<-f1st(y,x,p0=p0,kmn=kmn,sub=sub,kex=kex,inr=inr,xinr=xinr,qq=qq)
				if(a[[1]][1,1]>0){
					li<-length(a[[1]][,1])
					ind1<-(1:li)[a[[1]][,1]>0]
					ind1<-a[[1]][ind1,1]
					li<-length(ind1)
					tmp<-double(kexmx)
					if(kexmx<li+1) stop("kexmx too small")
					tmp[1]<-sum(a[[2]]^2)
					tmp[2:(li+1)]<-sort(ind1)
					covch<-rbind(covch,tmp)
					for(j in 1:li){
						kex0<-rbind(kex0,c(kex,ind1[j]))
					}
			
				}	
			}
			if(max(kex0)==0){mm<-m+1}
			else{mm<-mm+1
				ik<-length(kex0[,1])
				kexx<-kex0[2:ik,]
			}
		}
		tmp<-apply(covch,2,max)
		ind<-(1:kexmx)[tmp>0]
		li<-length(ind)

		covch<-covch[,ind]
		covch<-unique(covch)
		nc<-li
		nr<-length(covch)/li	
		covch<-matrix(covch,ncol=li,nrow=nr)

		lro<-length(covch[,1])
		covch<-as.double(covch)
		lco<-length(covch)/lro
		covch<-matrix(covch,nrow=lro,ncol=lco)
		if(lro>=2){
			rnk<-rank(covch[,1])
			lro<-1:lro
			lro[rnk]<-lro
			covch<-covch[lro,]
		}		
		list(covch)
}
