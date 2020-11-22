#' Selects the subsets specified by fmch. It is called by fmch. All subsets which are a subset of a specified  subset are removed. The remaining subsets are ordered by the sum of squares of the residuals
#'
#' @param nv The subsets specified by fmch
#' @param k  The number of covariates
#' @return ind The selected subsets.
#' @examples 
#' nv<-c(650,1962,160,5,3,4,577,1839,734)
#' nv<-matrix(nv,ncol=3)
#' a<-fselect(nv,6)
fselect<-function(nv,k){
	n<-length(nv[,1])
	ind<-1:n
	for(i in 1:(n-1)){
		tmpi<-decode(nv[i,1],k)[[2]]
		j<-i+1
		while(j<=n){
			if(nv[j,2]==nv[i,2]){j<-j+1}
			else{tmpj<-decode(nv[j,1],k)[[2]]
				tmpij<-tmpi*tmpj
				if(sum(abs(tmpij-tmpi))==0){
					ind[i]<- -i
					j<-n
				}
				j<-j+1
			}
		}
	}
	ind<-(1:n)[ind>0]
	list(ind)
}
