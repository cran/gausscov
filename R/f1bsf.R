#' Calculates non-parametric approximation to data y using fgenbsf.
#'
#' @param n Sample size
#' @param k Length of smallest interval
#' @param lam Factor for increasing  size of intervals
#' @param pr Proportional incerease in size of sample to reduce edge effects
#' @param mm  Parameter of fgentrig for the number of trigonometric functions
#' @param p0 Gaussian P-value threshold
#' @return ff The approximation
#' @return res The residuals.
#' @examples
#' data(vardata)
#' a<-f1bsf(vardata[,70],4,1.1,pr=0,mm=10)
f1bsf<-function(y,k,lam,pr=0.5,mm=20,p0=0.01){
	n<-length(y)
	m<-round(pr*n)
	nn<-n+2*m
	yy<-double(nn)
	if(m==0){yy<-y}
	if(m>0){
		yy[1:m]<-y[1]
		yy[(m+n+1):nn]<-y[n]
		yy[(m+1):(m+n)]<-y
	}
	x<-fgenbsf(nn,k,lam)[[1]]
	if(mm>0){
		x1<-fgentrig(nn,mm)[[1]]
		x<-cbind(x,x1)
	}
	a<-f1st(yy,x,p0=p0,sub=F)
	yy<-yy[(m+1):(m+n)]
	res<-a[[2]][(m+1):(m+n)]
	ff<-yy-res
	list(ff,res)
} 