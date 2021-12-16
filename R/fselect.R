#' It is called by fasb.R and frasb.R. All subsets which are a subset of a specified  subset are removed. The remaining subsets are ordered by the sum of squares of the residuals (fasb) or the scale (frasb)
#'
#' @param nv The subsets specified by fasb.R or frasb.R
#' @param k  The number of covariates
#' @return ind The selected subsets.
#' @example
#' b<-fasb(redwine[,12],redwine[,1:5],sel=FALSE)[[1]]
#' a<-fselect(b,11)[[1]]
#' b[a,]
fselect<-function(nv,k){
	n<-length(nv[,1])
	ir<-rank(-nv[,2],ties.method="first")
	ind<-1:n
	ind[ir]<-ind
	ind<-1:n
	rmv<-integer(n)
	for(i in 1:(n-1)){
		ii<-ind[i]
		tmpi<-decode(nv[ii,1],k)[[2]]
		ii2<-nv[ii,2]
		j<-i+1
		while(j<=n){
			jj<-ind[j]
			if(rmv[jj]==0){
				jj2<-nv[jj,2]
				if(jj2<ii2){
					tmpj<-decode(nv[jj,1],k)[[2]]
					ik<-1
					sbst<-1
					while(ik<=k){
						if(tmpj[ik]>tmpi[ik]){sbst<-0
							ik<-k+1
						}
						ik<-ik+1
					}
					if(sbst==1){rmv[jj]<-1}
				}
			}
			j<-j+1
		}
	}
	ind<-(1:n)[rmv==0]
	list(ind)
}
