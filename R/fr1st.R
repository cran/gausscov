#' Robust stepwise selection of covariates 
#'
#' @param y         Dependent variable.
#' @param x         Covariates.
#' @param cn        Constant  for Huber's psi-function.
#' @param cnr      The constants for Hampel's three part redescending psi function.
#' @param p0       The P-value cut-off.
#' @param sg       Scale value of residuals.
#' @param kmx     Specified number of included covariates.
#' @param mx      The maximum number of included covariates if the option sub=TRUE is used.
#' @param kex     The excluded covariates.
#' @param sub     Logical If TRUE best subset selected. 
#' @param inr       Logical if TRUE to include intercept 
#' @param xinr     Logical if TRUE intercept already included 
#' @param red     Logical It true Hampel's three part redescending psi function
#' @return pv        In order the subset ind, the regression coefficients, the Gaussian P-values, the standard P-values.
#' @return res   The residuals
#' @return stpv  The stepwise regression results: covariate, P-value and scale. 
#' @examples 
#' data(boston)
#' a<-fr1st(boston[,14],boston[,1:13],kex=7:8)
fr1st<-function(y,x,cn=1,cnr=c(1,3,5),p0=0.01,sg=0,kmx=0,mx=21,kex=0,sub=T,inr=T,xinr=F,red=F){
	if(mad(y)==0){stop("MAD=0")}
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
	y<-matrix(y,ncol=1)
	k<-length(x[1,])
	kex<-matrix(kex,nrow=1)
	lkx<-length(kex)
	if(xinr){
		if(kex[1]==0){q<-k-1}
		else{q<-k-1-lkx}
	}
	if(!xinr){
		if(inr){
			tmpx<-double(n)+1
			x<-cbind(x,tmpx)
			x<-matrix(x,nrow=n)
			xinr<-TRUE
			inr<-FALSE
		}
		q<-k
	}

	if(sg==0){sg<-mad(y)}
	kk<-length(x[1,])
	kex<-matrix(kex,nrow=1)
	lkx<-length(kex)
	if(kex[1]==0){q<-kk
		kex<-integer(kk+1)
	}
	else{q<-kk-lkx
		tmpi<-integer(kk+1)
		tmpi[kex]<-1
		kex<-tmpi
	}
	if(xinr){q<-q-1}
	kmax<-min(n,kk)
	km1<-kmax+1
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
		as.double(p0),
		as.integer(kmax),
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
		double(km1),
		as.logical(red),
		as.double(cnr),
		as.integer(kmx)
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
		ss01<-tmp[[23]][1:kmax]
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
			indd<-ind[1:(kmax-1)]
			pv<-frpval(y,x,cn=cn,q=q,ind=ind,inr=F,xinr=xinr)[[1]]
		}
		else{
			pv<-frpval(y,x,cn=cn,q=q,ind=ind,inr=F,xinr=F)[[1]]}
			li<-length(ind)
			if(xinr){pv[li,1]<-0
		}
		if(sub&(kmax>=2)){
			if(kmax>mx){stop("kmax too large (> mx) for all subset search")}
			if(xinr){
                			kk<-length(ind)
				sbsts<-frasb(y,x,cn=cn,cnr=cnr,p0=p0,q=q,ind=ind,inr=F,xinr=T,red=red)[[1]]
			}
			else{
				sbsts<-frasb(y,x,cn=cn,cnr=cnr,p0=p0,q=q,ind=ind,inr=F,xinr=F,red=red)[[1]]
			}
			sbsts<-matrix(sbsts,ncol=3)
			if(sbsts[1,1]>1){
				kk<-length(ind)
				tmv<-decode(sbsts[1,1],kk)[[1]]
				ind<-ind[tmv]
				if(xinr){
                				pv<-frpval(y,x,cn=cn,cnr=cnr,q=q,ind=ind,inr=T,xinr=F,red=red)
					res<-pv[[2]]
					pv<-pv[[1]]
				}
				else{
                				pv<-frpval(y,x,cn=cn,cnr=cnr,q=q,ind=ind,inr=F,xinr=F,red=red)
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
	list(pv,sg,res,stpv)
}

