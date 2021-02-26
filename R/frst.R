#' Robust stepwise selection of covariates 
#'
#' @param y         Dependent variable.
#' @param x         Covariates.
#' @param cn        Parameter for Huber's psi-function.
#' @param cnr      The constants for Hampel's three part redescending psi function.
#' @param p0       The P-value cut-off.
#' @param sg       Scale value of residuals.
#' @param nu       The order for calculating the P-value.
#' @param km      The maximum number of included covariates.
#' @param mx      The maximum number of included covariates if the option sub=TRUE is used.
#' @param kx       The excluded covariates.
#' @param sub     Logical If TRUE best subset selected. 
#' @param inr       Logical if TRUE to include intercept 
#' @param xinr     Logical if TRUE intercept already included 
#' @param red     Logical It true Hampel's three part redescending psi function
#' @return pv        In order the subset ind, the regression coefficients, the Gaussian P-values, the standard P-values.
#' @return res   The residuals
#' @return stpv  The stepwise regression results: covariate, P-value and scale. 
#' @examples 
#' data(boston)
#' a<-frst(boston[,14],boston[,1:13],km=10,sub=T)
frst<-function(y,x,cn=1,cnr=c(1,2,4),p0=0.01,sg=0,nu=1,km=0,mx=21,kx=0,sub=F,inr=T,xinr=F,red=F){
	if(mad(y)==0){stop("MAD=0")}
	if(sg==0){sg<-mad(y)}
	cnn<-cn
	if(red){cnn<-cnr[1]}
	tmpx<-cnn*(1:1000)/1000
	cpp<-sum(tmpx^2*dnorm(tmpx))*cnn/1000+cnn**2*(1-pnorm(cnn))
	cpp<-2*cpp
    	n<-length(y)
	x<-matrix(x,nrow=n)
	y<-matrix(y,ncol=1)
	k<-length(x[1,])
	if((!xinr)&inr){tmpx<-double(n)+1
		x<-cbind(x,tmpx)
		x<-matrix(x,nrow=n)
		xinr<-TRUE
	}
	kk<-length(x[1,])
	if(length(kx)==1){
		if(kx==0){lkx<-0}
	}
	else{lkx<-length(kx)}
	q<-kk-lkx
	if(xinr){q<-q-1}
	kex<-integer(kk+1)
	if(lkx==1){
		if(kx>0){
			kex[1]<-kx
		}
	}
	if(lkx>1){kex[1:lkx]<-kx}
	if((km>0)&inr){km=km+1}
	p00<-2
	if(km==0){
		p00<-p0
		km<-min(n,kk)
	}
	km1<-km+1
	tmp<-.Fortran(
		"robstepwise",		
		as.double(y),
		as.double(x),
		as.integer(n),
		as.integer(kk),	
		double(n*km1),
		double(n),
		double(n),
		double(km1),
		integer(kk+1),
		as.double(p00),
		as.integer(km),
		double(2*km1),
		double(km1),
		as.double(cn),
		as.double(cpp),
		as.double(sg),
		double(n),
		double(n),
		double(n),
		as.integer(km1),
		as.integer(kex),
		as.logical(xinr),
		as.double(nu),
		double(km1),
		as.logical(red),
		as.double(cnr)
	)
	kmax<-tmp[[11]]
	if(kmax==0){
		pv<-matrix(c(-2,0,0,0),nrow=1)
		res<-0
		stpv<-0
	}
	else if((kmax==1)&xinr){
		pv<-matrix(c(-1,0,0,0),nrow=1)
		res<-0
		stpv<-0
	}
	else{
		ss01<-tmp[[24]][1:kmax]
		stpv<-tmp[[12]]
		stpv<-matrix(stpv,ncol=2)
		stpv<-stpv[1:kmax,]
		sg<-tmp[[16]]
       		res<-tmp[[17]]
		stpv<-cbind(stpv,ss01)
		stpv<-matrix(stpv,ncol=3)
		ind<-stpv[1:kmax,1]
		if(xinr){ints<-ind[1]
			ind[1:(kmax-1)]<-ind[2:kmax]
			ind[kmax]<-ints
		}
		pv<-frrgp(y,x[,ind],cn=cn,q=q,inr=F,xinr=xinr)
		pv<-pv[[1]]
		pv[,1]<-ind
		li<-length(ind)
		if(xinr){pv[li,1]<-0}
		if(sub&(kmax>=2)){
			if(kmax>mx){stop("kmax too large (> mx) for all subset search")}
			if(xinr){
                			kk<-length(ind)
                			ind<-ind[1:(kk-1)]
				sbsts<-frmch(y,x,cn=cn,cnr=cnr,p0=p0,q=q,ind=ind,inr=T,xinr=F,red=red)[[1]]
			}
			else{
				sbsts<-frmch(y,x,cn=cn,cnr=cnr,p0=p0,q=q,ind=ind,inr=F,xinr=F,red=red)[[1]]
			}
			if(sbsts[1,1]>1){
				tmv<-decode(sbsts[1,1],k)[[1]]
				ind<-ind[tmv]
				if(xinr){
					q<-q+1
                				pv<-frrgp(y,x,cn=cn,cnr=cnr,q=q,ind=ind,inr=T,xinr=F,red=red)
					res<-pv[[2]]
					pv<-pv[[1]]
				}
				else{
                				pv<-frrgp(y,x,cn=cn,cnr=cnr,q=q,ind=ind,inr=F,xinr=F,red=red)
					res<-pv[[2]]
					pv<-pv[[1]]
				}
			}
		}
	}
	if((length(res)==n)&(xinr)){
		li<-length(pv[,1])
		pv[li,1]<-0
		stpv[1,1]<-0
	}
	list(pv,res,stpv)
}

