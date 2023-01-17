#'  Calculation of dependence graph using Gaussian stepwise selection.
#'
#' @param x Matrix of covariates.
#' @param p0  Cut-off P-value. 
#' @param ind Restricts the dependent nodes to this subset
#' @param kmn The minimum number of selected covariates for each node irrespective of cut-off P-value.
#' @param kmx The maximum number of selected covariates for each node irrespective of cut-off P-value.
#' @param mx Maximum nunber of selected covariates for each node for all subset search. 
#' @param nedge The maximum number of edges.
#' @param inr Logical if TRUE include an intercept.
#' @param xinr Logical if TRUE intercept already included.
#' @return ned Number of edges
#' @return edg List of edges together with P-values for each  edge and proportional reduction of sum of squared residuals
#' data(boston)
#' a<-fgr1st(boston[,1:13],ind=3:6) 
fgr1st<-function(x,p0=0.01,ind=0,kmn=0,kmx=0,nedge=10^5,inr=T,xinr=F){
	n<-length(x[,1])
	k<-length(x)/n
	ind<-matrix(ind,nrow=1)
	if(ind[1]==0){
		if(xinr){ind<-1:(k-1)}
		else{ind<-1:k}
	}
	li<-length(ind)
	p0<-p0/k
	x<-matrix(x,nrow=n)
	inrr<-inr
	if((!xinr)&(inr)){tmpx<-double(n)+1
		x<-cbind(x,tmpx)
		k<-k+1
		xinr<-T
		inr<-F
	}
	xx<-x
	if(kmx==0){kmx<-min(n,k)}
	xinrr<-0
	if(xinr){xinrr<-1}
	tmp<-.Fortran(
		"graphst",
		as.double(xx),
		as.double(x),
		as.integer(n),
		as.integer(k),
		double(n),
		double(n),
		double(n),
		integer(k),
		as.double(p0),
		as.integer(kmx),	
		double((k+1)*2),
		integer(nedge*2),
		integer(1),
		integer(k),
		as.integer(xinrr),
		double(k),
		as.integer(nedge),
		double(k),
		as.integer(kmn),
		as.integer(li),
		as.integer(ind),
		double(nedge)
		)
	ned<-tmp[[13]]
	if(ned>0){
        		edg<-tmp[[12]]
		edg<-matrix(edg,nrow=nedge,ncol=2)
		edgp<-tmp[[22]]
		edg<-cbind(edg,edgp)
        		edg<-edg[1:ned,]
		edg<-matrix(edg,nrow=ned,ncol=3)
		eedg<-c(0,0,0)
		kk<-k
		if(xinr){kk<-k-1}
		for(i in 1:kk){
			tmpi<-(1:ned)[edg[,1]==i]
			ind<-edg[tmpi,2]
			li<-length(ind)
			if(li>0){
				qq<-k-li
				if(xinr){qq<-qq-1}
				ind1<-ind
				a<-fpval(x[,i],x,ind1,inr=inr,xinr=xinr)
				la<-length(a[[1]][,1])
				ind<-(1:la)[a[[1]][,1]>0]
				ind<-(a[[1]][ind,1])
				li<-length(ind)
				for(ii in 1:li){
					eedg0<-cbind(i,ind[ii],a[[1]][ii,3])
					eedg<-rbind(eedg,eedg0)
				}	
			}
		}
	}
    	else{
       		edg<-integer(3)
		eedg<-matrix(edg,nrow=1,ncol=3)
	}
	edg<-eedg[2:(ned+1),]
	list(ned,edg)
}
