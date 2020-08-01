#' Selects the subsets specified by fmch. It is called by fmch. All subsets which are a subset of a specified  subset are removed. The remaining subsets are ordered by the sum of squares of the residuals
#'
#' @param nv The subsets specified by fmch
#' @param k  The number of covariates
#' @return ind The selected subsets.
#' @examples 
#' nv<-c(650,1962,160,1033,394,1730,577,1839,334,895)
#' a<-fselect(nv,11)
fselect<-function(nv,k){
	n<-length(nv)
	ind<-1:n
	for(i in n:2){
		tmpi<-decode(nv[i],k)[[2]]
		j<-i-1
		while(j>=1){
			tmpj<-decode(nv[j],k)[[2]]
			tmpij<-tmpi*tmpj
			if(sum(abs(tmpij-tmpi))==0){
				ind[i]<- -i
				j<-0
			}
			j<-j-1
		}
	}
	ind<-(1:n)[ind>0]
	list(ind)
}
