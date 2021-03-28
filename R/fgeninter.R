#' Generates interactions of a given order from the covariates
#'
#' @param x Covariates
#' @param ord     Order of interactions
#' @param inr Logical to include intercept
#' @return xx  All interactions of order at most ord.
#' @examples 
#' data(boston)
#' bosint<-fgeninter(boston[,1:13],2)[[1]]
fgeninter<-function(x,ord,inr=T){
	n<-length(x[,1])
	if(inr){
		tmpx<-double(n)+1
		x<-cbind(tmpx,x)
	}
	k<-length(x[1,])
	kk<-choose(k-1+ord,k-1)
	tmp<-.Fortran(
		"genint",
		as.double(x),
		double(n*kk),
		as.integer(n),
		as.integer(k),
		as.integer(kk),
		as.integer(ord),
		integer(ord)
		)
	xx<-tmp[[2]]
	dim(xx)<-c(n,kk)
	xx<-xx[,2:kk]
	tmp<-double(n)+1
	xx<-cbind(xx,tmp)
	xx<-matrix(xx,ncol=kk)
	list(xx)
}


