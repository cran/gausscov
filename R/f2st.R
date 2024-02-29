#' Repeated stepwise disjoint selection of covariates 
#'
#' @param y Dependent variable.
#' @param x Covariates.
#' @param p0  The cut-off P-value.
#' @param kmn The minimum number of included covariates irrespective of cut-off P-value.
#' @param kmx The maximum number of included covariates irrespective of cut-off P-value.
#' @param kex The excluded covariates.
#' @param mx  The maximum number covariates for an all subset search.
#' @param lm  The maximum number of linear approximations. 
#' @param sub Logical, if T select best subset.
#' @param inr Logical if T include an intercept .
#' @param xinr Logical if T intercept already included.
#' @param qq   The number of covariates to choose from. If qq=-1 the number of covariates of x is used.
#' @return cpv   In order the linear approximation,  the included covariates, the Gaussian P-values.
#' @examples 
#' data(boston)
#' bostint<-fgeninter(boston[,1:13],2)[[1]]
#' a<-f2st(boston[,14],bostint,lm=3,sub=FALSE)
f2st<-function(y,x,p0=0.01,kmn=0,kmx=0,kex=0,mx=21,lm=9^9,sub=T,inr=T,xinr=F,qq=-1){
        if(xinr&(kmn==1)){stop("only intersect left")}
        n<-length(y)
        k<-length(x)/n
        x<-matrix(x,nrow=n,ncol=k)
        y<-matrix(y,ncol=1)
        lkx<-length(kex)
        kex<-matrix(kex,nrow=1)
        if(!xinr){
                if(inr){
                        tmpx<-double(n)+1
                        x<-cbind(x,tmpx)
                        k<-k+1
                        x<-matrix(x,nrow=n,ncol=k)
                        xinr<-TRUE
                        inr<-FALSE
                }
        }
        if(xinr){
                if(lkx==1){
                        if(kex[1]==0){kex[1]<-k}
                }
                kex<-c(kex,k)
        }
        if(xinr){
                if(kmn>0){kmn<-kmn+1}
                if(kmx>0){kmx<-kmx+1}
        }
        lkx<-length(kex)
        qq<-k-lkx
        pv<- matrix(c(0,0,0,0,0),nrow=1)
        kv<-1
        mn<-0
        while(kv>0.5){
                a<-f1st(y,x,p0=p0,kmn=kmn,kmx=kmx,kex=kex,mx=mx,sub=sub,inr=inr,xinr=xinr,qq=qq)[[1]]
                if(a[1,1]<=0){kv<-0}
                if(kv>0.5){
                        lv<-length(a[,1])
                        llv<-(1:lv)[a[,1]>0]
                        lllv<-length(llv)
                        cv<-a[1:lllv,1]
                        pv<-a[1:lllv,3]
                        mn<-mn+1
                        mnt<-integer(lllv)+mn
                        mnt<-matrix(mnt,ncol=1)
                        if(mn==1){
                                cpv<-cbind(mnt,cv,pv)
                        }
                        else{
                                cpvn<-cbind(mnt,cv,pv)
                                cpv<-rbind(cpv,cpvn)
                        }
                        kex<-c(kex,cv)
                        kv<-1
                        if(mn>=lm){kv<-0}
                        if(length(kex)>=k-1){kv<-0}
                }
        }
        list(cpv)
}
