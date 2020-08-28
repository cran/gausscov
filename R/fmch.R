#' Calculates all possible subsets and selects those where each included covariate is significant. 
#'
#' 
#'
#' @param y The dependent variable.
#' @param x  The covariates.
#' @param p0 Cut-off p-value for significance.
#' @param q The number of covarites from which to choose 
#' @param ind   The indices of subset of covariates for which all subsets are to be considered.
#' @param select Logical. If TRUE remove all subsets of chosen sets.
#' @param inr Logical If TRUE include intercept.
#' @param xinr Logical If TRUE intercept already included.

#' @return nvv Coded list of subsets with number of covariates and sum of squared residuals
#' @examples 
#' data(redwine)
#' a<-fmch(redwine[,12],redwine[,1:11])
fmch<-function(y,x,p0=0.01,q=-1,ind=0,sel=T,inr=T,xinr=F){
	n<-length(y)
	x<-matrix(x,nrow=n)
	k<-length(x)/n
	    ind<-matrix(ind,nrow=1)
	if((!xinr)&inr){
		tmpx<-double(n)+1
		x<-cbind(x,tmpx)
		x<-matrix(x,nrow=n)
		xinr<-TRUE
		k<-k+1
        if(ind[1]>0){ind<-c(ind,k)}
	}
	if(ind[1]>0){
		x<-x[,ind]
		x<-matrix(x,nrow=n)
		k<-length(x)/n
	}
	else{ind<-1:k}
	ss<-double(2^k+2)
	q<-max(q,k)
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
		double(k^2),
		integer(k),
		as.logical(xinr),
		double(2^k+2),
		integer(2^(k+2)),
		as.double(ss),
		as.double(p0),
		as.integer(q)
	)
	ss<-tmp[[17]]
	inv<-(1:(2^k+2))[ss>0]
	ss<-ss[inv]
	nv<-tmp[[16]]
	dim(nv)<-c(2^(k+1),2)
	nv<-matrix(nv,ncol=2)
	nv<-nv[inv,]
	llv<-length(inv)
	if(llv>0){
		nv<-tmp[[16]]
		nv<-matrix(nv,ncol=2)
		ind<-rank(ss,ties.method="first")
		inv<-1:llv
		inv[ind]<-inv
		ss<-ss[inv]
		nv<-nv[inv,]
		nv<-matrix(nv,ncol=2)
		nv1<-nv[,1]
		nv2<-nv[,2]
#
#	select approximations 
#
		nvv<-cbind(nv1,nv2,ss)
		nvv<-matrix(nvv,ncol=3)
		if(sel&(llv>1)){
			nv2<-nvv[,2]
			ind<-rank(nv2,ties.method="first")	
			inv<-1:llv		
			inv[ind]<-inv
			indd<-fselect(nvv,k)[[1]]
			nv1<-nv1[indd]
			nv2<-nv2[indd]
			ss<-ss[indd]
			nvv<-cbind(nv1,nv2,ss)
			nvv<-matrix(nvv,ncol=3)
		}
	}
	else{nvv<-matrix(c(-1,0),nrow=1)}
	list(nvv)
}

