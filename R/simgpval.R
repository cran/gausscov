#' Simulates Gaussian P-values
#'
#' @param y The dependent variable
#' @param x The chosen subset of covariates
#' @param i The chosen covariate
#' @param nsim The number of simulations
#' @param qq The total number of covariates available. If qq=-1 the number of covariates of x is taken.
#' @param plt   Logical if TRUE plot Gaussian P-values
#' @return pgv  The Gaussian P-value of the ith covariate and the relative frequency with which the Gaussian covariates are better. 
#' @examples 
#' data(snspt)
#'snspt<-matrix(snspt,nrow==3253,ncol=1)
#' a<-flag(snspt,3253,1,12)
#' simgpval(a[[1]],a[[2]],7,10,plt=FALSE)
simgpval<-function(y,x,i,nsim,qq=-1,plt=TRUE){
        n<-length(y)      
        k<-length(x)/n
        if(qq==-1){qq<-k}    
        a<-lm(y~x) 
        ssx<-sum(a$res^2)   
        px<-summary(a)[[4]][(i+1),4] 
        px<-1-(1-px)^(qq-k+1)
        res<-double(nsim) 
        if(k>1){xi<-x[,-i]}
        ic<-0     
        for(isim in 1:nsim){
                RSS<-sum(y^2)
                res[isim]<-1
                for(j in 1:(qq-k+1)){
                        Z<-rnorm(n) 
                        if(k>1){
                                b<-lm(y~I(xi)+I(Z)) 
                        }
                        else{
                                b<-lm(y~I(Z))
                        } 
                        if(j==1){res[isim]<-summary(b)[[4]][(k+1),4]}      
                        RSS<-pmin(RSS,sum(b$res^2))
                }      
                 if(RSS<ssx){ic<-ic+1} 
        }
        if(plt){plot(sort(res),t="l")}  
        rfx<-ic/nsim               
        pgv<-c(px,rfx)                 
        list(pgv)
}


