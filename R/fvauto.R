#' Vector autoregressive approximation
#'
#' @param x  Variables
#' @param n  Sample size
#' @omx      Maximum lag
#' @param p0 The cut-off P-value.
#' @param kmn Minimum number of included covariates irrespective of cut-off P-value
#' @param kmx Maxmum number of included covariates irrespective of cut-off P-value
#' @param mx  The maximum number covariates for an all subset search.
#' @param kex The excluded covariates.
#' @param sub Logical, if TRUE best subset selected.
#' @param inr Logical, if TRUE include intercept if not present.
#' @return res The selected lagged variables for each variable 
#' @return res2 The regression coefficients and P-values
#' @return res4 The residuals
#' @examples 
#' data(abcq)
#' a<-fvauto(abcq,240,10)
fvauto<-function(x,n,omx,p0=0.01,kmn=0,kmx=0,mx=21,kex=0,sub=T,inr=TRUE){
	k<-length(x)/n
	res<-integer((k+1)*5)
	res0<-res
	res3<-res
	res2<-list(res)
	res4<-list(1:n)
	for(i in 1:k){
		xl<-flag(x,n,i,omx)
		a<-f1st(xl[[1]],xl[[2]],p0=p0,kmn=kmn,kmx=kmx,mx=mx,kex=kex,sub=sub,inr=inr)
		if(a[[1]][1,1]>0){
			b<-a[[2]]
			res4<-c(res4,list(b))
			ss1<-sum(b^2)
			a<-a[[1]]
			la<-length(a[,1])
			res1<-res0
			ind1<-(1:la)[a[,1] >0]
			lind1<-length(ind1)
			res1[1]<-i
			res1[2:(lind1+1)]<-a[ind1,1]
			res<-rbind(res,res1)
			ss<-double(4)
			ss[1]<-ss1
			a<-rbind(a,ss)

			res2<-c(res2,list(a))
#			res4<-c(res4,list(b))
		}
		else{
			res<-rbind(res,res0)
			res2<-c(res2,list(res0))
			res4<-c(res4,list(res0))
		}
	}
	nr<-length(res)/((k+1)*5)
	res<-matrix(res,nrow=nr,ncol=(k+1)*5)
	j<-1
	while(j<=(k+1)*5){
		tmp<-sum(res[,j])
		if(tmp>0){jj<-j
			j<-j+1
		}
		else{jj<- j-1
			j<-(k+1)*5+1
		}
	}
	res<-res[,1:jj]
	res<-res[2:nr,]
	lres2<-length(res2)
	res2<-res2[2:lres2]
	res4<-matrix(res4,nrow=k+1)
	res4<-res4[2:(k+1),]
	list(res,res2,res4)
}
			
			
		
	
