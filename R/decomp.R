#' Decompose a given coded interaction ic into its component parts
#'
#' @param ic The numbers of the interactions.
#' @param k  The number of covariates of x in fgeninter.R with no intercept, number plus 1 if inr=TRUE in fgeninter.R
#' @param ord The order of the interactions.
#' @param inc The indices of the interaction covariates with no dummy covariates when all powers are calculated including dummy covariates. This is returned by fgeninter.R.
#' @return decom The component parts of the interactions.
#' @examples 
#' bosint<-fgeninter(boston[,1:13],3,4)
#' a<-decomp(100,14,3,inc=bosint[[2]])
decomp<-function(ic,k,ord,inc=0){
    	ic<-ic+1
	m<-length(ic)
	inc<-matrix(inc,nrow=1)
	if(inc[1]>0){ic<-inc[ic]}
	tmp<-.Fortran(
		"degenint",
		as.integer(ic),
		as.integer(m),
		as.integer(k),
		as.integer(ord),
		integer(m*ord),
		integer(ord)
		)
	dim(tmp[[5]])<-c(m,ord)
	decom<-tmp[[5]]
	list(decom)
}
