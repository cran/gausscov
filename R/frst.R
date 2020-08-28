#' Robust stepwise selection of covariates 
#'
#' @param y         Dependent variable
#' @param x         Covariates
#' @param cn        Parameter for Huber's psi-function
#' @param p0       The P-value cut-off
#' @param sg       Scale value of residuals
#' @param nu       The order for calculating the P-value
#' @param km      The maximum number of included covariates
#' @param mx      The maximum number of included covariates if the option subset =TRUE is used.
#' @param kx       The excluded covariates
#' @param sb     Logical If TRUE best subset selected 
#' @param of Logical if TRUE to include intercept 
#' @param xof Logical if TRUE intercept already included 
#' @return pv    In order the subset ind, the regression coefficients, the P-values, the standard P-values.
#' @return res   The residuals
#' @return stpv  The stepwise regression results: covariate, P-value and scale. 
#' @examples 
#' data(boston)
#' a<-frst(boston[,14],boston[,1:13],km=10,sub=T)
frst<-function(y,x,cn=1,p0=0.01,sg=0,nu=1,km=0,mx=20,kx=0,sub=F,inr=T,xinr=F){
	if(mad(y)==0){stop("MAD=0")}
	if(sg==0){sg<-mad(y)}
	tmpx<-cn*(1:1000)/1000
	cpp<-sum(tmpx^2*dnorm(tmpx))*cn/1000+cn**2*(1-pnorm(cn))
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
	q<-kk
	kex<-integer(kk+1)
	lke<-length(kx)
	if(lke==1){
		if(kx>0){
			kex[1]<-kx
		}
	}
	if(lke>1){kex[1:lke]<-kx}
	if((km>0)&inr){km=km+1}
	p00<-2
	if(km==0){
		p00<-p0
		km<-min(n,kk)
	}
	km1<-km+1
	pv<-double(3*kk)
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
		double(km1)
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
		pv<-frobregp(y,x[,ind],cn=cn,q=q,inr=F,xinr=xinr)
		pv<-pv[[1]]
		pv[,1]<-ind
		li<-length(ind)
		if(xinr){pv[li,1]<-0}
		if(sub&(kmax>=2)){
			if(kmax>mx){stop("kmax too large (> mx) for all subset search")}
			if(xinr){
                		kk<-length(ind)
                		ind<-ind[1:(kk-1)]
				sbsts<-frmch(y,x,cn=cn,p0=p0,ind=ind,inr=T,xinr=F)[[1]]
			}
			else{
				sbsts<-frmch(y,x,cn=cn,p0=p0,ind=ind,inr=F,xinr=F)[[1]]
			}
			if(sbsts[1,1]>1){
				tmv<-decode(sbsts[1,1],k)[[1]]
				ind<-ind[tmv]
				if(xinr){
                			pv<-frobregp(y,x,cn=cn,q=q,ind=ind,inr=T,xinr=F)
					res<-pv[[2]]
					pv<-pv[[1]]
				}
				else{
                			pv<-frobregp(y,x,cn=cn,q=q,ind=ind,inr=F,xinr=F)
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

