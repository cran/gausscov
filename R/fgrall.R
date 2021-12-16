#'  Calculation of dependency graph using Gaussian all subset selection.
#'
#' @param x Matrix of covariates.
#' @param p0  The cut-off P-value.
#' @param kmx The maximum number of included covariates for each node irrespective of cut-off P-value.
#' @param mx The maximum number of covariates.
#' @param inr Logical, if TRUE  to include intercept.
#' @param xinr Logical, if TRUE  intercept already included.
#' @return ned Number of edges
#' @return edg List of edges with Gaussian P-value and percentage of sum of squared residuals explained by edge
#' a<-fgrall(redwine[,1:8])
fgrall<-function(x,p0=0.01,kmx=0,mx=21,inr=T,xinr=F){
	n<-length(x[,1])
	k<-length(x)/n
	if(k>mx) stop("number of covariates too large k > mx")
	p0<-p0/k
	x<-matrix(x,nrow=n)
	if(inr){tmpx<-double(n)+1
		x<-cbind(x,tmpx)
		k<-k+1
		xinr<-T
		inr<-F
	}
	qq<-k-1
	kk<-k
	if(xinr){qq<-k-2
		kk<-k-1
	}
	edge<-c(0,0,0,0)
	for(i in 1:kk){
		tmp<-integer(k)
		tmp[i]<-1
		ind1<-(1:k)[tmp==0]
		li<-length(ind1)
		xx<-x[,ind1]
		kx<-length(xx)/n
		if(xinr){kx<-kx-1}
		a<-fasb(x[,i],xx,p0=p0,q=qq,inr=inr,xinr=xinr)[[1]]
		a<-matrix(a,ncol=3)
		if(a[1,1]>0){
			ind2<-decode(a[1,1],li)[[1]]
			li<-length(ind2)
			qqq<-kx-li+1
			pv<-fpval(x[,i],xx,ind=ind2,inr=F,xinr=xinr,q=qqq)
			pv<-pv[[1]]
			ll<-length(pv[,1])
			ind3<-ind1[pv[1:(ll-1),1]]
			pv[1:(ll-1),1]<-ind3
			if(ll>1){
				for(j in 1:(ll-1)){
					jj<-pv[j,1]
					pjj<-pv[j,3]
					prc<-pv[j,5]
					edge<-rbind(edge,c(i,jj,pjj,prc))
				}
			}
		}
	}
	ne<-length(edge[,1])
	if(ne>1){edge<-edge[2:ne,]}
	ne<-ne-1
	list(ne,edge)
}

