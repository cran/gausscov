#' Calculates all possible subsets where each included covariate is significant. Decode with fdecode.
#'
#' 
#'
#' @param y The dependent variable.
#' @param x  The covariates.
#' @param p0 Cut-off p-value for significance.
#' @param q The number of covariates from which to choose. Equals number of covariates minus length of ind if q=-1. 
#' @param ind   The indices of a subset of covariates for which all subsets are to be considered.
#' @param sel  Logical. If TRUE calls fselect.R to remove all subsets of chosen sets.
#' @param inr Logical If TRUE include intercept.
#' @param xinr Logical If TRUE intercept already included.

#' @return nvv Coded list of subsets with number of covariates and sum of squared residual
#' @examples 
#' data(redwine)
#' a<-fasb(redwine[,12],redwine[,1:11])
fasb<-function(y,x,p0=0.01,q=-1,ind=0,sel=T,inr=T,xinr=F){
	n<-length(y)
	x<-matrix(x,nrow=n)
	k<-length(x)/n
	if(q==-1){
		q<-k
		if(xinr){q<-k-1}
	}
	ind<-matrix(ind,nrow=1)
	if(ind[1]==0){
		if(!xinr){
			if(inr){
				tmpx<-double(n)+1
				x<-cbind(x,tmpx)
				x<-matrix(x,nrow=n)
				k<-k+1
				xinr<-TRUE
				inr<-FALSE
			}
		}
	}
	else{
		if(xinr){ind<-c(ind,k)}
		x<-x[,ind]
		if(inr){
			tmpx<-double(n)+1
			x<-cbind(x,tmpx)
			x<-matrix(x,nrow=n)
			xinr<-TRUE
			inr<-FALSE
		}
	}
	k<-length(x)/n
	
	kmxx<-2^k
	if(xinr){kmxx<-2^(k-1)}
	tmp<-.Fortran(
		"lmmdch",
		as.double(y),
		as.double(x),
		as.integer(n),
		as.integer(k),	
		double(n*k),
		double(n*k),
		double(n),
		double(n),
		double(k),
		double(k),
		double(k),
		double(k**2),
		integer(k),
		as.logical(xinr),
		double(2**k+2),
		integer(2**(k+1)),
		double(2**k+2),
		as.double(p0),
		as.integer(q)
	)
	ss<-tmp[[17]]
	inv<-(1:(2**k+2))[ss>0]
	ss<-ss[inv]
	rs<-rank(ss)
	nv<-tmp[[16]]
	dim(nv)<-c(2**k,2)
	nv<-matrix(nv,ncol=2)
	nv<-nv[inv,]
	nv<-matrix(nv,ncol=2)
	llv<-length(inv)
	if(llv>0){
		ind1<-rank(ss,ties.method="first")
		inv<-1:llv
		inv[ind1]<-inv
		ss<-ss[inv]
		nv<-nv[inv,]
		nv<-matrix(nv,ncol=2)
		nv1<-nv[,1]
		nv2<-nv[,2]
#
#	select approximations 
#
		nvv<-cbind(nv1,nv2,ss)
		nr<-length(nvv)/3
		nvv<-matrix(nvv,ncol=3,nrow=nr)
		if(sel&(llv>1)){
			ind2<-rank(-nv2,ties.method="first")	
			inv<-1:llv		
			inv[ind2]<-inv
			nvv<-nvv[inv,]
			kk<-k
			if(xinr){kk<-k-1}
			indd<-fselect(nvv,kk)[[1]]
			li<-length(indd)
			if(li==1){nvv<-nvv[li,]}
			else{li<-1:li
				nvv<-nvv[indd,]
				rl<-rank(nvv[,3])
				li[rl]<-li
				nvv<-nvv[li,]
			}
		}
		if(llv==1){nvv<-cbind(nv1[1],nv2[1],ss[1])
			nvv<-matrix(nvv,ncol=3,nrow=1)
		}
	}
	else{nvv<-matrix(c(0,0,0),ncol=3,nrow=1)}
	nr<-length(nvv)/3
	nvv<-matrix(nvv,ncol=3,nrow=nr)
	list(nvv)
}

