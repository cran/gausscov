#' Decodes the number of a subset selected by fasb.R to give the  covariates 
#'
#'
#' @param ns The number of the subset.
#' @param k The number of covariates excluding the intercept.
#' @return ind The list of covariates.
#' @return set A binary vector giving the covariates.
#' @examples 
#' decode(19,8)
decode<-function(ns,k){
	set<-integer(k)
	tmp<-.Fortran(
		"decode",
		as.integer(ns),
		as.integer(k),
		integer(k)
		)
	set<-tmp[[3]]
	ind<-(1:k)[set==1]
	ind<-matrix(ind,nrow=1)
	list(ind,set)
}
