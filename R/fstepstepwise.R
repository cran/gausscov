#' Repeated stepwise selection of covariates 
#'
#' @param y Dependent variable
#' @param x Covariates
#' @param alpha  The P-value cut-off
#' @param nu The order statistic of Gaussian covariates used for comparison
#' @param kmax The maximum number of included covariates
#' @param kexk The excluded covariates
#' @param lmax  The maximum number of linear approximations 
#' @param intercept Logical to include intercept 
#' @param chkintercept Logical to include or exclude intercept dependent on the P-value
#' @param misclass Logical giving the number of misclassifications if appropriate,  eg for binary y 
#' @return pv In order, the number of linear approximation, the included covariates, the P-values, sum of squared residuals and if appropriate number of misclassifications.
#' @examples 
#' data(colonx)
#' data(colony)
#' fstepstepwise(colon.y,colon.x,0.01,lmax=10,misclass=TRUE) 
fstepstepwise<-function(y,x,alpha,nu=1,kmax=0,kexk=0,lmax=10^10,intercept=TRUE,chkintercept=FALSE,misclass=FALSE){
	xx<-x
	kk<-length(x[1,])+1
	n<-length(y)
	k<-length(x[1,])
	indn<-double(k)+1
	ln<-length(indn)
	pv<-double(4)
	dim(pv)<-c(1,4)
	if(misclass){
		pv<-double(5)
		dim(pv)<-c(1,5)
	}
	ic<-0
	kv<-1
	kex<-0
	nm<-0
	is<-0
	while(kv>0.5){
		tmp<-fstepwise(y,x,alpha,nu,kmax,kexk,intercept=intercept,chkintercept=chkintercept,misclass=misclass)[[1]]
		lv<-length(tmp)
		kv<-tmp[1,1]
		if(kv>0){
			if(lv>3){lv<-length(tmp[,1])}
			else{lv<-1}
			nm<-nm+1
			kexk<-c(kexk,tmp[1:lv,1])
			ic<-ic+1
			kx<-double(4*lv)
			dim(kx)<-c(lv,4)
			if(misclass){
				kx<-double(5*lv)
				dim(kx)<-c(lv,5)
			}
			kx[1:lv,1]<-ic
			kx[1:lv,2]<-tmp[1:lv,1]
			kx[1:lv,3]<-tmp[1:lv,2]
			kx[1:lv,4]<-tmp[1:lv,3]
			if(misclass){
				kx[1:lv,5]<- tmp[1:lv,4]
			}
			pv<-rbind(pv,kx)
			is<-is+1
		}
		if(nm==lmax){kv<-0}
	}
	lk<-length(pv)/4
	if(misclass){
		lk<-length(pv)/5
	}		
	if(lk>1){pv<-pv[2:lk,]
		if(misclass){
		dim(pv)<-c(lk-1,5)
		}
		else{dim(pv)<-c(lk-1,4)}
		ind<-(1:(lk-1))[pv[1:(lk-1),2]>0]
		pv<-pv[ind,]
	}
	list(pv)
}
