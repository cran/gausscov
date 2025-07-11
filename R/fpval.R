#' Calculates the regression coefficients, the P-values and the standard P-values for the chosen subset ind
#
#' @param y The dependent variable.
#' @param x  The covariates.
#' @param ind The indices of the subset of the covariates whose P-values are required.
#' @param inr Logical If TRUE intercept to be included
#' @param xinr Logical If TRUE intercept already included.
#' @param qq   The total  number of covariates from which ind was selected. If qq=-1 the number of covariates of x minus length ind plus 1 is taken.
#' @return apv In order the subset ind, the regression coefficients, the Gaussian P-values, the standard F P-values and the proportion of sum of squares explained.
#' @return res The residuals.
#' @examples 
#' a<-fpval(boston[,14],boston[,1:13],c(1,2,4:6,8:13))
fpval<-function(y,x,ind,inr=T,xinr=F,qq=-1){
        n<-length(y)
        kx<-length(x)/n
        li<-length(ind)
        ki<-length(ind)
        x<-matrix(x,nrow=n)
        y<-matrix(y,ncol=1)
        ind<-matrix(ind,nrow=1)
        if(!xinr){
                if(inr){
                        tmpx<-double(n)+1
                        x<-cbind(x,tmpx)
                        xinr<-TRUE
                }
        }
        kx<-length(x)/n
        if(xinr){if(max(ind)<kx){ind<-c(ind,kx)}}
        kii<-length(ind)
        if(xinr){ind<-sort(ind)}
        xx<-x[,ind]
        xx<-matrix(xx,nrow=n)
        b<-lm(y~0+xx)
        ind1<-(1:kii)[is.na(b$coef)==FALSE]
        ind<-ind[ind1]
        res<-as.double(b$res)
        result<-summary(b)[[4]]
        kii<-length(result[,1])
        apv<-double(4*kii)
        apv<-matrix(apv,nrow=kii,ncol=4)
#       tvl<-matrix(result[1:kii,3],ncol=1)
#	print(dim(apv))
#	print(length(ind))
        apv[,1]<-ind
        apv[,2]<-matrix(result[1:kii,1],ncol=1)
        apv[,4]<-matrix(result[1:kii,4],ncol=1)
        apv[,3]<-apv[,4]
#       tvl<-1-1/(1+tvl^2/(n-ki))
        if(qq==-1){qq<-kx-li+1
		if(xinr){qq<-qq-1}
        }
        apv[1:kii,3]<-pbeta(apv[1:kii,4],1,qq)
        if(xinr){
	apv[kii,3]<-apv[kii,4]
                apv[kii,1]<-0
        }
        list(apv,res)
}
