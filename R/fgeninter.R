#' Generates interactions of a given order from the covariates
#'
#' @param x Covariates
#' @param ord   Order of interactions
#' @return xx  All interactions of order at most ord excluding 0-1 covariates
#' @return intx Decomposes a given interaction covariate of xx
#' @examples 
#' data(boston)
#' bosint<-fgeninter(boston[,1:13],3)
fgeninter<-function(x,ord){
	n<-length(x[,1])
	tmpx<-double(n)+1
	x<-cbind(x,tmpx)
	k<-length(x[1,])
	kk<-choose(k-1+ord,k-1)
	intx<-integer(kk*ord)
	tmp<-.Fortran(
		"genint",
		as.double(x),
		double(n*kk),
		as.integer(n),
		as.integer(k),
		as.integer(kk),
		as.integer(intx),
		as.integer(ord),
		integer(ord),
		integer(1)
		)
	xx<-tmp[[2]]
	xx<-matrix(xx,ncol=kk,nrow=n)
	kkx<-tmp[[9]]
	xx<-xx[,1:kkx]
	intx<-tmp[[6]]
	intx<-matrix(intx,ncol=ord,nrow=kk)	
	list(xx,intx)
}


