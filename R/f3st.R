#' Stepwise selection of covariates after exclusion of already selected covariates 
#'
#' @param y Dependent variable.
#' @param x Covariates.
#' @param m The number of iterations
#' @param p0 The cut-off P-value.
#' @param kmn The minimum number of included covariates irrespective of cut-off P-value.
#' @param kmx The maximum number of included covariates irrespective of cut-off P-value.
#' @param kex The excluded covariates.
#' @param mx  The maximum number covariates for an all subset search.
#' @param sub Logical, if TRUE best subset selected.
#' @param inr Logical, if TRUE include intercept if not present.
#' @param xinr Logical, if TRUE intercept already included 
#' @qq   The number of covariates to choose from. If qq=-1 the number of covariates of x is used.
#' @param kexmx The maximum number of covariates in an approximation
#' @return covch The standard deviation of residuals and the selected covariates ordered in increasing size of the standatrd deviations.
#' @return lai The number of rows of covch
#' @examples 
#' data(leukemia)
#' a<-f3st(leukemia[[1]],leukemia[[2]],m=1,kmn=5,kexmx=10)
f3st<-function(y,x,m,mxa,p0=0.01,kmn=0,kmx=0,kex=0,mx=21,sub=T,inr=T,xinr=F,qq=-1,kexmx=100){
        covch<-double(kexmx)
        covch<-matrix(covch,nrow=1,ncol=kexmx)
        a<-f1st(y,x,p0=p0,kmn=kmn,kmx=kmx,mx=mx,kex=kex,sub=sub,inr=inr,xinr=xinr,qq=qq)
        if(a[[1]][1,1]>0){
                ali<-length(a[[1]][,1])
                ind<-(1:ali)[a[[1]][,1]>0]
                ind<-a[[1]][ind,1]
                ali<-length(ind)
                covch[1,1]<-sd(a[[2]])
                a<-a[[1]]
                ind<-sort(ind)
                covch[1,2:(ali+1)]<-ind
                li<-length(ind)
		kex<-ind
                lai<-1
                if(m>0){
                        a<-f3sti(y,x,covch,ind,m,mxa,kexmx=kexmx,p0=p0,kmn=kmn,kmx=kmx,mx=mx,kex=kex,sub=sub,inr=inr,xinr=xinr,qq=qq)
                        covch<-a[[1]]
                        lai<-length(covch[,1])
                }
                if(lai>0){
                        tmp<-apply(covch,2,max)
                        len<-length(tmp)
                        ind<-(1:len)[tmp>0]
                        covch<-covch[,ind]
                }
        }
        else{
                covch<-0
                lai<-0
        }
        covch<-round(10^10*covch)*10^-10
        covch<-unique(covch)
        list(covch,lai)
}
