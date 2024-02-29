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
#' @param qq   The number of covariates to choose from. If qq=-1 the number of covariates of x is used.
#' @return ned Number of edges
#' @return edg List of edges together with P-values for each  edge and proportional reduction of sum of squared residuals
#' data(boston)
#' a<-fgr1st(boston[,1:13],ind=3:6) 
fgr1st<-function(x,p0=0.01,ind=0,kmn=0,kmx=0,mx=21,nedge=10^5,inr=T,xinr=F,qq=-1){
        dx<-dim(x)
        n<-dx[[1]]
        m<-dx[[2]]
        edg<-c(0,0,0)
        if(min(ind)==0){ind<-1:m}
        x<-x[,ind]
        li<-length(ind)
        for(i in 1:li){
                gr<-f1st(x[,i],x,p0=p0,kmn=kmn,kmx=kmx,kex=i,mx=mx,sub=T,inr=inr,xinr=xinr,qq=qq)
                if(gr[[1]][1,1] >=1){
                        lgr<-length(gr[[1]][,1])
                        il<-(1:lgr)[gr[[1]][,1] >0]
                        knt<-integer(length(il))+i
                        gr1<-gr[[1]][il,1]
                        pv1<-gr[[1]][il,3]
                        edg1<-cbind(knt,gr1,pv1)
                        edg<-rbind(edg,edg1)
                }
        }
        ne<-length(edg)/3
        edg<-edg[2:ne,]
        ne<-ne-1
	for(i in 1:ne){
		edg[i,1]<-ind[edg[i,1]]
		edg[i,2]<-ind[edg[i,2]]
	}
        list(ne,edg)
}
