#' Vector autoregressive approximation
#'
#' @param x  Variables
#' @n	     Sample size
#' @omx      Maximum lag
#' @param p0 The cut-off P-value.
#' @return res The selected lagged variables for each variable 
#' @examples 
#' data(abcq)
#' a<-fvauto(abcq,240,10)
fvauto<-function(x,n,omx,p0=0.01){
	k<-length(x)/n
	xl<-flag(x,n,omx)[[2]]
	p0<-p0/k
	res<-integer((k+1)*5)
	res0<-res
	for(i in 1:k){
		x1<-x[(omx+1):n,i]
		a<-f1st(x1,xl)[[1]]
		la<-length(a[,1])
		res1<-res0
		if(la>=2){
			ind1<-a[1:(la-1),1]

			res1[1]<-i
			res1[2:la]<-ind1
		}
		else{
			res1[1]<-i
		}
		res<-rbind(res,res1)
	}
	nr<-length(res)/((k+1)*5)
	res<-matrix(res,nrow=nr,ncol=(k+1)*5)
	j<-1
	while(j<=(k+1)*5){
		tmp<-sum(res[,j])
		if(tmp>0){j<-j+1}
		else{jj<- j-1
			j<-(k+1)*5+1
		}
	}
	res<-res[,1:jj]
	res<-res[2:nr,]
	list(res)
}
			
			
		
	
