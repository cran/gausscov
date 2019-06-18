#' Simulates the number of false positives for given dimensions (n,k) and given order statistic nu
#'
#' @param n  The dimension of dependent variable
#' @param k  The number of covariates
#' @param alpha Cut-off p-value
#' @param nu Order statistic
#' @param kmax Maximum number of false positives
#' @param nsim Number of simulations
#' @return res Histogram of number of false positives.
#' @return mn Mean number of false positives.
#' @return ss Standard deviation of number of false positives
#' @examples 
#' fsimords(100,1000,0.05,3,6,nsim=100)
fsimords<-function(n,k,alpha,nu,kmax,nsim=100){
	set.seed(11051944)
	y<-double(n)
	y[1]<-1
	ss<-0
	res<-double(2*kmax)
	dim(res)<-c(kmax,2)
	res[1:kmax,1]<-0:(kmax-1)
	mn<-0
	ss<-0
	for(isim in 1:nsim){
		tmpx<-rnorm(n*k)
		dim(tmpx)<-c(n,k)
		tmp<-fstepwise(y,tmpx,alpha,kmax,nu=nu)[[1]]
		if(tmp[1,1]==0){res[1,2]<-res[1,2]+1}
		if(tmp[1,1]>0){
			kl<-length(tmp[,1])
			mn<-mn+kl
			ss<-ss+kl**2
			if(kl<=kmax-1){
				res[kl+1,2]<-res[kl+1,2]+1
			}
			else{
				res[kmax,2]<-res[kmax,2]+1
			}
		}
	}
	mn<-mn/nsim
	ss<-ss/nsim
	ss<-sqrt(ss-mn**2)	
	res[,2]<-res[,2]/nsim
	ss<-sqrt(sum((res[,1]-mn)**2*res[,2]))
	list(res,mn,ss)
}
