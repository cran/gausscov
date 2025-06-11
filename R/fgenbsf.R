# Generates intervals with a basis function
#'
#' @param n Sample size
#' @param k Length of smallest interval
#' @param lam   Factor for increasing sizes of intervals
#' @return x The matrix of the intervals.
#' @examples 
#' a<-fgenbsf(100,4,1.2)
fgenbsf<-function(n,k,lam){
	j0<-ceiling(log(n/k)/lam)
	ii<-0
	for(j in 0:j0){
	      kk<-ceiling(lam^j*k)
	      i<-1
	      while(i*kk<=n){
	      	   a1<-(i-1)*kk+1
		   a2<-a1+kk-1
		   ii<-ii+1
               	   if(ceiling(a2+k/2)<=n){
			ii<-ii+1
		   }
		   i<-i+1
	      }   
	}
	ind<-matrix(nrow=ii,ncol=2)
	jj<-0
	for(j in 0:j0){
	      kk<-ceiling(lam^j*k)
	      i<-1
	      while(i*kk<=n){
	      	   a1<-(i-1)*kk+1
		   a2<-a1+kk-1
		   jj<-jj+1
		   ind[jj,1]<-a1
     	   	   ind[jj,2]<-a2
               	   if(ceiling(a2+k/2)<=n){
			jj<-jj+1
			ind[jj,1]<-floor(a1+k/2)
			ind[jj,2]<-ceiling(a2+k/2)
		   }
		   i<-i+1
	      }   
	}
	x<-matrix(nrow=n,ncol=jj)
	for(k in 1:jj){
	      n1<-ind[k,1]
	      n2<-ind[k,2]
	      li<-n2-n1+1
              tmpx<-double(li)
              tmpy<-(0:(li-1))/(li-1)
              for(ji in 1:li){
		       tmpx[ji]<-tmpy[ji]^3*(1-tmpy[[ji]])^3
	      }
	      x[1:n,k]<-0
	      x[n1:n2,k]<-tmpx
	 }
	list(x)
}

		