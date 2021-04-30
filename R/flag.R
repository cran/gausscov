#' Calculation of lagged covariates 
#'
#' @param x Covariates.
#' @param n The sample size.
#' @param lag The maximum lag.
#' @return  y The first covariate without a lag, the dependent covariate.
#' @return xl The lagged covariates with lags of order 1:lag starting with the first covariate.
#' @examples 
#' data(abcq)
#' abcql<-flag(abcq,240,16)
#' a<-f1st(abcql[[1]],abcql[[2]])
flag<-function(x,n,lag){
	k<-length(x)/n
	x<-matrix(x,nrow=n,ncol=k)
	tmp<-.Fortran(
		"lagg",
		as.double(x),
		as.integer(n),
		as.integer(k),
		as.integer(lag),
		double((n-lag)*k*lag),
		double(n-lag)
		)
	xl<-tmp[[5]]
	xl<-matrix(xl,nrow=n-lag,ncol=k*lag)
	y<-tmp[[6]]
	list(y,xl)
}
