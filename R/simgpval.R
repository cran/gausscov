#' Simulates Gaussian P-values
#'
#' @param y The dependent variable 
#' @param x The covariates
#' @param i The chosen covariate
#' @param nsim The number of simulations
#' @param plt   Logical if TRUE plot Gaussian P-values
#' @return pgv  Th P-value of the ith covariates and the relative frequency with which the Gaussian covariates are better. 
#' @examples 
#' data(snspt)
#' a<-flag(snspt,3253,12)
#' simgpval(a[[1]],a[[2]],7,10,plt=FALSE)
simgpval<-function(y,x,i,nsim,plt=TRUE){
	n<-length(y)                    # number of observations
	q<-length(x[1,])                # number of covariates
	a<-lm(y~x) 		# regress y on x
	
	ssx<-sum(a$res^2)              # sum of squared residuals
	px<-summary(a)[[4]][(i+1),4]    # standard F P-value of ith covariate of x
	res<-double(nsim)               # P-values of the Gaussian covariates
	xi<-x[,-i]                      # x without the ith covariate x_i
	ic<-0                           # counts how often the Gaussian covariates  are better than x_i
	for(isim in 1:nsim){
      		Z<-rnorm(n)                         # a Gaussian covariate
      		b<-lm(y~I(xi)+I(Z))                 # x_i replaced by Z
      		res[isim]<-summary(b)[[4]][(q+1),4] # p-value of Z
      		RSS<-sum(b$res^2)                   # sum of squared residuals with x_i   replaced by Z
     		 if(RSS<ssx){ic<-ic+1}               # if Z better than x_i  increase   ic by 1
	}
	if(plt){plot(sort(res),t="l")}      # plot ordered Gaussian P-values
	rfx<-ic/nsim               
 	pgv<-c(px,rfx)                 		# print P-value of x_i and relative frequency 
	list(pgv)
}


