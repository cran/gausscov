#Stepwise selection of covariates based on Gaussian P-values
#'
#' @param y Dependent variable.
#' @param x Covariates.
#' @param p0 The cut-off P-value.
#' @param kmn The minimum number of included covariates irrespective of cut-off P-value.
#' @param kmx The maximum number of included covariates irrespective of cut-off P-value.
#' @param kex The excluded covariates.
#' @param mx  The maximum number covariates for an all subset search.
#' @param sub Logical, if TRUE best subset selected.
#' @param inr Logical, if TRUE include intercept if not present.
#' @param xinr Logical, if TRUE intercept already included 
#' @param qq   The number of covariates to choose from. If qq=-1 the number of covariates of x is used.
#' @return pv In order the included covariates, the regression coefficient values, the Gaussian P-values, the standard P-values. If pv(1)=-1, no subset returned.
#' @return res The residuals.
#' @examples 
#' data(boston)
#' bostint<-fgeninter(boston[,1:13],2)[[1]]
#' a<-f1st(boston[,14],bostint,kmn=10,sub=TRUE)
f1st<-function(y,x,p0=0.01,kmn=0,kmx=0,kex=0,mx=21,sub=T,inr=T,xinr=F,qq=-1){
        if(xinr&(kmn==1)){stop("only intersect left")}
	n<-length(y)
	if(length(x)/n==1){x<-matrix(x,nrow=n,ncol=1)}
        dm<-dim(x)
        n<-dm[1]
        k<-dm[2]
        lkx<-length(kex)
        kex<-matrix(kex,nrow=1)
        if(!xinr){
                if(inr){
                        tmpx<-double(n)+1
                        x<-cbind(x,tmpx)
                        k<-k+1
                        x<-matrix(x,nrow=n,ncol=k)
                        xinr<-TRUE
                        inr<-FALSE
                }
        }
        if(xinr){
                if(kmn>0){kmn<-kmn+1}
                if(kmx>0){kmx<-kmx+1}
        }
	a<-integer(k)
	if(xinr){a[k]<-1}
	if(qq==-1){qq<-k}
	if(min(kmx)==0){kmx<-k}
        lkx<-length(kex)
        tmp<-.Fortran(
                "fstepwise",
                as.double(y),
                as.double(x),
                as.integer(n),
                as.integer(k),
                double(n),
                double(n),
                as.integer(a),
                as.double(p0),
                as.integer(kmx),
                double(2*(k+1)),
                as.integer(kex),
                double(n),
                double(k),
                as.integer(qq),
                as.integer(kmn),
                as.integer(lkx)
        )
	kmax<-tmp[[9]]
        if(kmax==0){
                pv<-matrix(c(-1,0,0,0),nrow=1)
                res<-0
                stpv<-0
        }
        if(kmax==1){
		if(xinr){
                	pv<-matrix(c(-1,0,0,0),nrow=1)
                	res<-0
                	stpv<-0
		}
		else{
			stpv<-tmp[[10]]
               			stpv<-matrix(stpv,ncol=2)
			stpv<-stpv[1,]
			stpv<-matrix(stpv,ncol=2,nrow=1)
			i<-stpv[1,1]
			b<-fpval(y,x,ind=i,inr=inr,xinr,qq=qq)
			pv<-b[[1]]
			res<-b[[2]]
		}
	}
        if(kmax>=2){
                stpv<-tmp[[10]]
                stpv<-matrix(stpv,ncol=2)
                stpv<-stpv[1:(kmax+1),]
                stpv<-matrix(stpv,ncol=2)
                ind<-stpv[1:kmax,1]
		ind<-sort(ind)
               	li<-length(ind)
		b<-fpval(y,x,ind=ind,inr=inr,xinr=xinr,qq=qq)
		pv<-b[[1]]
		res<-b[[2]]
		k1<-length(b[[1]][,1])
		if(xinr){k1<-k1-1
			pmx<-max(b[[1]][1:k1,3])
			if(pmx<p0){sub=FALSE}
		}
		if((mx<k1)&(sub==TRUE)){sub<-FALSE}
                if(sub==FALSE){
                   	b<-fpval(y,x,ind=ind,inr=inr,xinr=xinr,qq=qq)
                   	pv<-b[[1]]
	   		res<-b[[2]]
       	 	}
                if(sub==TRUE){
                 	if(sub&(li>=2)){
                                	sbsts<-fasb(y,x,p0=p0,ind=ind,inr=inr,xinr=xinr,qq=qq)[[1]]
                        		if(sbsts[1,1]>0){
                                		nss<-sbsts[1,1]
                                		tmv<-decode(nss,li)[[1]]
                                		ind<-ind[tmv]
                                		if(xinr){if(max(ind)<k){ind<-c(ind,k)}}
                                		a<-fpval(y,x,ind,inr=inr,xinr=xinr,qq=qq)
                                		li<-length(ind)
                                		res<-a[[2]]
                                		pv<-a[[1]]
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
        }
	 list(pv,res)
}

