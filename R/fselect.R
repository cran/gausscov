#' Selects the subsets specified by flmmdch. It is called by flmmdch
#'
#' All subsets which are a subset of a specified  subset are removed. The remaining subsets are ordered by the sum of squares of the residuals
#'
#' @param nv The subsets specified by flmmdch
#' @param k  The variables
#' @return ind The selected subsets.
#' @examples 
#' fselect(nv,13)
fselect<-function(nv,k){
	n<-length(nv[,1])
	ind<-1:n
	jj<-1
	tmp<-decode(nv[jj,1],k)
	for(i in ind[1:(n-1)]){
		if(i>0){tmpi<-decode(nv[i,1],k)
			for(j in ind[(i+1):n]){
				if(j>0){tmpj<-decode(nv[j,1],k)
					tmpij<-tmpi*tmpj
					if(sum(abs(tmpij-tmpj))==0){ind[j]<- -j}
				}
			}
		}
	}
	ind<-ind[ind>0]
	list(ind)
}
