#' Simulates the number of false positives for given dimensions (n,k) and given order statistics nu
#'
#' @param n  The dimension of dependent variable.
#' @param k  The number of covariates.
#' @param p0 Cut-off p-value.
#' @param nu Order statistics.
#' @param km Maximum number of selected covariates.
#' @param nsim Number of simulations.
#' @return p   Histogram of number of false positives.
#' @return mn Mean number of false positives.
#' @return ss Standard deviation of number of false positives.
#' @examples 
#' a<-fsimords(100,100,0.01,c(5,10),15,nsim=100)
fsimords<-function(n,k,p0,nu,km,nsim=500){
	y<-double(n)
	y[1]<-1
	y[2]<--1
	lnu<-length(nu)
	res<-double(lnu*nsim)
	res<-matrix(res,nrow=lnu)
	for(isim in 1:nsim){
		tmpx<-rnorm(n*k)
		dim(tmpx)<-c(n,k)
		tmp<-f1st(y,tmpx,km=km,inr=FALSE)[[3]]
		pv<-tmp[,2]
		for(j  in 1:lnu){
			res[j,isim]<-0
			i<-1
			while(i <=km){
				pp<-qbeta(1-pv[i],k+2-i,1)
				pp<-1-pbeta(pp,k+3-i-nu[j],nu[j])
				if(pp<p0){
					res[j,isim]<-res[j,isim]+1
				}
				else{i<-km}
				i<-i+1
			}
		}
	}
	p<-double(lnu*(km+1))
	p<-matrix(p,nrow=lnu)
	mn<-double(lnu)
	ss<-double(lnu)
	ind<-0:km
	for(j in 1:lnu){
		for( i in 0:km){
			p[j,i+1]<-length((1:nsim)[res[j,]==i])/nsim
		}
		mn[j]<-sum(ind*p[j,1:(km+1)])
		ss[j]<-sqrt(sum(ind^2*p[j,1:(km+1)])-mn[j]^2)
	}
	list(p,mn,ss)
}
