#' Calculates all possible valid approximations which are of subsets of the given subset.
#'
#'
#'
#' @param y The dependent variable.
#' @param x  The covariates.
#' @param p0 Cut-off p-value for significance.
 #' @param ind   The indices of the subset of covariates for which all subsets are to be considered.
#' @param inr Logical If TRUE include intercept.
#' @param xinr Logical If TRUE intercept already included.
#' @param qq The number of covariates from which to choose. Equals number of covariates minus length of ind if qq=-1.
#' @return nvv Coded list of subsets with number of covariates and sum of squared residual
#' @examples 
#' data(redwine)
#' a<-fasb(redwine[,12],redwine[,1:11])
fasb<-function(y,x,p0=0.01,ind=0,inr=T,xinr=F,qq=-1){
          dx<-dim(x)
          n<-dx[1]
          k<-dx[2]
          ind<-matrix(ind,nrow=1)
          if(ind[1]==0){ind<-1:k}
          if(!xinr){
		if(inr){
	     	tmpx<-double(n)+1
	     	x<-cbind(x,tmpx)
	     	k<-k+1
	     	ind<-c(ind,k)
	     	inr<-F
	     	xinr<-T
	     	}
          }
          if(qq==-1){qq<-k
		if(xinr){qq<-qq-1}
          }  
          x<-x[,ind]
          k<-length(x)/n
          kmxx<-2^k
#          if(xinr){kmxx<-2^(k-1)}
          xinrr<-0
          if(xinr){xinrr<-1}
          tmp<-.Fortran(
	                "lmmdch",
	                as.double(y),
	                as.double(x),
	                as.integer(n),
	                as.integer(k),
	                double(n*k),
	                double(n*k),
	                double(n),
	                double(n),
	                double(k),
	                double(k),
	                double(k),
	                double(k**2),
	                integer(k),
	                as.integer(xinrr),
	                double(2**k+2),
	                integer(2**(k+1)),
	                double(2**k+2),
	                as.double(p0),
	                as.integer(qq)
          )
          ss<-tmp[[17]]
          inv<-(1:(2**k+2))[ss>0]
          ss<-ss[inv]
          rs<-rank(ss)
          nv<-tmp[[16]]
          dim(nv)<-c(2**k,2)
          nv<-matrix(nv,ncol=2)
          nv<-nv[inv,]
          nv<-matrix(nv,ncol=2)
          llv<-length(inv)
          if(llv>0){
                ind1<-rank(ss,ties.method="first")
                inv<-1:llv
                inv[ind1]<-inv
                ss<-ss[inv]
                nv<-nv[inv,]
                nv<-matrix(nv,ncol=2)
                nv1<-nv[,1]
                nv2<-nv[,2]
#
#       select approximations 
#
                nvv<-cbind(nv1,nv2,ss)
                nr<-length(nvv)/3
                nvv<-matrix(nvv,ncol=3,nrow=nr)
                if(llv==1){nvv<-cbind(nv1[1],nv2[1],ss[1])
                        nvv<-matrix(nvv,ncol=3,nrow=1)
                }
          }
          else{nvv<-matrix(c(0,0,0),ncol=3,nrow=1)}
          nr<-length(nvv)/3
          nvv<-matrix(nvv,ncol=3,nrow=nr)
          list(nvv)
}

