#' Decodes the number of a subset selected by flmmdch to give the  covariates 
#'
#'
#' @param j The number of the subset
#' @param k The number of covariates
#' @return set A binary vector giving the covariates
#' @examples 
#' decode(19,8)
decode<-function(j,k){
	jj<-j
	set<-double(k)
	for(i in (k-1):0){
		if(2^i<=jj){
			set[i+1]<-1
			jj<-jj-2^i
		}
	}
	return(set)
}
