#' Stepwise selection of covariates 
#'
#' @param y Dependent variable
#' @param x Covariates
#' @param alpha  The P-value cut-off
#' @param nu The order statistic of Gaussian covariates used for comparison
#' @param kmax The maximum number of included covariates
#' @param kexk The excluded covariates
#' @param intercept Logical to include intercept 
#' @param chkintercept Logical to include or exclude intercept dependent on the P-value
#' @param misclass Logical The number of misclassifications if appropriate,  eg for binary y 
#' @return pv The selected covariates in order together with P-values, sum of squared residuals and if appropriate number of misclassifications.
#' @examples 
#' data(colonx)
#' data(colony)
#' fstepwise(colon.y,colon.x,0.01)
#' fstepwise(colon.y,colon.x,1,kmax=4) 
fstepwise<-function(y,x,alpha,nu=1,kmax=0,kexk=0,intercept=TRUE,chkintercept=FALSE,misclass=FALSE){
	n<-length(y)
	x<-as.matrix(x)
	y<-as.matrix(y)
	k<-length(x[1,])
	yy<-y
	if(kmax==0){kmax<-n}
	if(intercept){tmpx<-double(n)+1
		x<-cbind(tmpx,x)
	}
	k<-length(x[1,])
	lke<-length(kexk)
	kex<-double(k+1)
	kex[1:lke]<-kexk
	kmax1<-kmax+1
	res<-y
	pv<-double(3*k)
	dim(pv)<-c(k,3)
	if(misclass){
		pv<-double(4*k)
		dim(pv)<-c(k,4)
	}
	tmp<-.Fortran(
		"fststepwise",
		as.double(yy),
		as.double(x),
		as.integer(n),
		as.integer(k),	
		double(n),
		as.double(res),
		integer(k+1),
		as.double(alpha),
		as.integer(kmax),
		double(2*kmax1),
		as.integer(kmax1),
		as.integer(kex),
		as.logical(intercept),
		as.double(nu),
		double(kmax1),
		as.logical(chkintercept),
		as.logical(misclass),
		integer(kmax1)
	)
	kkmax<-tmp[[9]]
	pvv<-tmp[[10]]
	dim(pvv)<-c(kmax1,2)
	if(kkmax>0){
		minss<-tmp[[15]]
		is<-1
		kkkmax<-kkmax
		if(pvv[1,1]==0){is<-2
			kkkmax<-kkmax-1
		}
		pv[1:kkkmax,1]<-pvv[is:kkmax,1]
		pv[1:kkkmax,2]<-pvv[is:kkmax,2]
		pv[1:kkkmax,3]<-minss[is:kkmax]
		if(misclass){
			pv[1:kkkmax,4]<-tmp[[18]][is:kkmax]
		}
		pv<-pv[1:kkkmax,]
		res<-tmp[[6]]
		if(misclass){ 
			pv<-as.matrix(pv,nrow=kkkmax,ncol=4)
		}	
		else{
			pv<-as.matrix(pv,nrow=kkkmax,ncol=3)
		}
		if(kkmax==1){pv<-t(pv)}
	}
	else{pv<-c(0,0)
		dim(pv)<-c(1,2)
	}
	list(pv)
}
