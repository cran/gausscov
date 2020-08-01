#' Calculates all possible subsets and selects those where each included covariate is significant using a robustified version of fmch.R . 
#'
#' @param y The dependent variable.
#' @param x  The covariates.
#' @param cn The constant for Huber's psi function.
#' @param p0 Cut-off p-value.
#' @param sg The scale parameter.
#' @param ind The subset for which the results are required.
#' @param sl Logical. If TRUE remove all subsets of chosen sets.
#' @param of Logical If TRUE to inlude intercept. 
#' @param xof Logical If TRUE x already has intercept.
#' @return nv List of subsets with number of covariates and scale.
#' @examples 
#' data(boston)
#' a<-frmch(boston[,14],boston[,1:8]) 
frmch<-function(y,x,cn=1,p0=0.01,sg=0,ind=0,sl=T,of=T,xof=F){
	if(mad(y)==0){stop('mad(y) is zero')}
	if(sg==0){sg<-mad(y)}
	tmpx<-cn*(1:1000)/1000
	cpp<-sum(tmpx^2*dnorm(tmpx))*cn/1000+cn**2*(1-pnorm(cn))
	cpp<-2*cpp
	n<-length(y)
	x<-matrix(x,nrow=n)
	k<-length(x)/n
	if(!xof){q<-k}
	else{q<-k-1}
	if((!xof)&of){
		tmpx<-double(n)+1
		x<-cbind(x,tmpx)
		x<-matrix(x,nrow=n)
		xof<-TRUE
		k<-k+1
	}
	ind<-matrix(ind,nrow=1)
        if(ind[1]>0){
                x<-x[,ind]
                x<-matrix(x,nrow=n)
                k<-length(x)/n
        }
        else{ind<-1:k}

	tmp<-.Fortran(
		"roblmmdch",
		as.double(y),
		as.double(x),
		as.integer(n),
		as.integer(k),
		as.double(p0),	
		double(n*k),
		double(n*k),
		double(n*k),
		double(n),
		double(k),
		double(k),
		double(k),
		double(k^2),
		double(n),
		double(k),
		as.double(cn),
		as.double(sg),
		double(3),
		as.double(cpp),
		integer(k),
		integer(k),
		double(2*k),
		as.logical(xof),
		integer(2^(k+2)),
		double(2^k),
		as.integer(q)
	)
	ss<-tmp[[25]]
	inv<-(1:2^k)[ss>0]
	ss<-ss[inv]
	llv<-length(inv)
	if(length(inv)==0){nvv=c(-1,-1,-1)
			nvv<-matrix(nvv,nrow=1)
	}
	if(length(inv)==1){
			tmp<-decode(inv,k)[[2]]
			nv<-sum(tmp)
			nvv<-c(inv,nv,ss[inv])
			nvv<-matrix(nvv,nrow=1)
	}
	if(length(inv) >1){
		nv<-tmp[[24]]
		dim(nv)<-c(2^(k+1),2)
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
		if(sl&(llv>1)){
			indd<-fselect(nv1,k)[[1]]
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


