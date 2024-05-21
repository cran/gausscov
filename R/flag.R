#' Calculation of lagged covariates 
#' 
#' @param x The matrix of covariates.
#' @praam n The sample size
#' @praam i The selected covariate
#' @param lag The maximum lag.

#' @return  y The ith  covariate of x  without a lag, the dependent covariate.
#' @return xl The lagged covariates with lags of order 1:lag starting with the first covariate.
#' @examples 
#' data(abcq)
#' abcql<-flag(abcq,240,1,16)
#' a<-f1st(abcql[[1]],abcql[[2]])
flag<-function(x,n,i,lag){
	m<-length(x)/n
	x<-matrix(x,nrow=n,ncol=m)
	tmp<-.Fortran(
		"lagg",
		as.double(x),
		as.integer(n),
		as.integer(m),
		as.integer(i),
		as.integer(lag),
		double((n-lag)*m*lag),
		double(n-lag)
		)
	y<-tmp[[7]]
	xl<-tmp[[6]]
	xl<-matrix(xl,nrow=n-lag,ncol=m*lag)
	list(y,xl)
}
