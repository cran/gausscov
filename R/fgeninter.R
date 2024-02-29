#' Generates interactions of a given order from the covariates
#'
#' @param x Covariates
#' @param ord   Order of interactions
#' @return xx  All interactions of order at most ord 
#' @return intx Decomposes a given interaction covariate of xx
#' @examples 
#' data(boston)
#' bosint<-fgeninter(boston[,1:13],2)
fgeninter<-function(x,ord){
        dx<-dim(x)
        n<-dx[1]
        k<-dx[2]
        tmpx<-double(n)+1
        x<-cbind(x,tmpx)
        dx<-dim(x)
        n<-dx[1]
        k<-dx[2]
        kk<-choose(k-1+ord,k-1)
        intx<-integer(kk*ord)
        xx<-double(n*kk)+0
        xx<-as.matrix(xx,nrow=n,ncol=kk)
        ji=0
        tmp<-.Fortran(
                "genint",
                as.double(x),
                as.double(xx),
                as.integer(n),
                as.integer(k),
                as.integer(kk),
                as.integer(intx),
                as.integer(ord),
                integer(ord+1),
                as.integer(ji)
                )
        xx<-tmp[[2]]
        xx<-matrix(xx,ncol=kk,nrow=n)
        kkx<-tmp[[9]]
        xx<-xx[,1:kkx]
        intx<-tmp[[6]]
        intx<-matrix(intx,ncol=ord,nrow=kk)
        list(xx,intx)
}


