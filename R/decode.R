#' Decodes the number of a subset selected by flmmdch to give the  covariates 
#'
#'
#' @param j The number of the subset
#' @param k The number of covariates including the intercept or larger
#' @return ind The list of covariates
#' @return set A binary vector giving the covariates
#' @examples 
#' decode(19,8)
decode<-function(j,k){
	jj<-j
	set<-double(k)
	for(i in (k-1):0){
		if(jj>=2^i){
			set[i+1]<-1
			jj<-jj-2^i
		}
	}
	ind<-(1:k)[set==1]
	ind<-matrix(ind,nrow=1)
	list(ind,set)
}
