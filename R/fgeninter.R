#' Generates interactions of a given order from the covariates
#'
#' @param x Covariates
#' @param ord   Order of interactions
#' @param inr Logical to include intercept
#' @param idv List of 0-1 dummy covariates
#' @return xx  All interactions of order at most ord excluding powers of 0-1 covariates
#' @return inc The indices of the covariates in xx when all powers are calculated including dummy variables to be used in decomp.R.
#' @examples 
#' data(boston)
#' bosint<-fgeninter(boston[,1:13],3,idv=4)
fgeninter<-function(x,ord,inr=TRUE,idv=0){
	n<-length(x[,1])
	tmpx<-double(n)+1
	x<-cbind(tmpx,x)
	k<-length(x[1,])
	kk<-choose(k-1+ord,k-1)
	dv<-integer(k)+0
	if(length(idv)==1){
		if(idv>0){dv[idv]<-1}
	}
	if(length(idv)>1){dv[idv]<-1}
	kex<-integer(kk)
	tmp<-.Fortran(
		"genint",
		as.double(x),
		double(n*kk),
		as.integer(n),
		as.integer(k),
		as.integer(kk),
		as.integer(kex),
		as.integer(ord),
		as.integer(dv),
		integer(ord),
		integer(k)
		)
	xx<-tmp[[2]]
	dim(xx)<-c(n,kk)
	xx<-xx[,2:kk]
	tmpx<-double(n)+1
	xx<-cbind(xx,tmpx)
	xx<-matrix(xx,ncol=kk)
	kex<-tmp[[6]]
	if(max(idv)>0){kex<-(1:kk)[kex==1]
		tmpi<-integer(kk)
		tmpi[kex]<-1
		tmpi<-(1:kk)[tmpi==0]
		xx<-xx[,tmpi]
		kk<-length(tmpi)
		inc<-tmpi
	}
	else{inc<-0}
	if(!inr){xx<-xx[,1:(kk-1)]}
	list(xx,inc)
}


