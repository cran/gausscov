#' Decompose a given coded interaction ic into its component parts
#'
#' @param ic The number of the interaction.
#' @param k  The number of covariates of x in fgeninter with no intercept.
#' @param ord The order of the interactions.
#' @return decom The component parts of the interaction.
#' @examples 
#' a<-decomp(7783,14,8)
decomp<-function(ic,k,ord){
    	ic<-ic+1
	m<-length(ic)
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
