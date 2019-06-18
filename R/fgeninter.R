#' generation of interactions 
#'
#' @param x Covariates
#' @param ord     Order of interactions
#' @param intercept Logical to include intercept
#' @return xx  All interactions of order at most ord.
#' @examples 
# data(gausscov)
#' bostinter<-fgeninter(boston[,1:13],7)[[1]]
#'dim(bostinter)
fgeninter<-function(x,ord,intercept=TRUE){
	n<-length(x[,1])
	if(intercept){
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
	list(xx)
}
