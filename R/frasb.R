#' Calculates all possible subsets and selects those where each included covariate is significant using a robustified version of fasb.R based on Huber's psi-function or Hampel's redescending three  psi-function . 
#'
#' @param y The dependent variable.
#' @param x  The covariates.
#' @param cn The constant for Huber's psi function.
#' @param cnr The constants for Hampel's three part redescending psi-function
#' @param p0 Cut-off p-value.
#' @param q  If q>0 the number of covariates from which ind was chosen.
#' @param sg The scale parameter.
#' @param ind The subset for which the results are required.
#' @param sel Logical. If TRUE removes all subsets of chosen sets.
#' @param inr  Logical If TRUE inlude intercept. 
#' @param xinr Logical If TRUE intercept included in x.
#' @param red  Logical If true Hampel's three part redescending psi function
#' @return nvv Coded list of subsets with number of covariates and scale ordered according to scale.
#' @examples 
#' data(boston)
#' a<-frasb(boston[,14],boston[,1:8]) 
frasb<-function(y,x,cn=1,cnr=c(1,3,5),p0=0.01,q=-1,sg=0,ind=0,sel=T,inr=T,xinr=F,red=F){
	if(mad(y)==0){stop('mad(y) is zero')}
	cnn<-cn
	if(red){cnn<-cnr[3]
		tmpx<-cnn*(1:1000)/1000
		cpp<-0
		for(i in 1:1000){
			cpp<-cpp+fpsired(tmpx[i],cnr)^2*dnorm(tmpx[i])
		}
		cpp<-cpp*cnn/1000
	}
	else{tmpx<-cnn*(1:1000)/1000
		cpp<-sum(tmpx^2*dnorm(tmpx))*cnn/1000+cnn**2*(1-pnorm(cnn))
	}
	cpp<-2*cpp
	n<-length(y)
	x<-matrix(x,nrow=n)
	k<-length(x)/n
	ind<-matrix(ind,nrow=1)
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
			}
		}
	}
	else{
		x<-x[,ind]
		if(!xinr){
			if(inr){
				tmpx<-double(n)+1
				x<-cbind(x,tmpx)
				x<-matrix(x,nrow=n)
				xinr<-TRUE
			}
		}
	}
	k<-length(x)/n
	if(sg==0){sg<-mad(y)}
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
		as.logical(xinr),
		integer(2^(k+2)),
		double(2^k),
		as.integer(q),
		double(n),
		as.logical(red),
		as.double(cnr)
	)
	ss<-tmp[[25]]
	inv<-(1:2^k)[ss>0]
	llv<-length(inv)
	if(llv==0){	nvv<-matrix(c(0,0,0),nrow=1)
	}
	else{ss<-ss[inv]
		rs<-rank(ss)
		nv<-tmp[[24]]
		nv<-matrix(nv,ncol=2)
		nv<-nv[inv,]
	}
	if(llv==1){
		nvv<-c(nv,ss[inv])
		nvv<-matrix(nvv,nrow=1)
	}
	if(llv>1){
#		ss<-ss[inv]
		nv1<-nv[,1]
		nv2<-nv[,2]
#
#	select approximations 
#
		nvv<-cbind(nv1,nv2,ss)
		nvv<-matrix(nvv,ncol=3)
		if(sel){
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
	}
	list(nvv)
}


