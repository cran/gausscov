#' Converts directed graph (a,b) not equal to (b,a) into an undirected graph (a,b)=(b,a)
#'
#' @param gr A directed graph
#' @return gr The undirected graph
#' @examples 
#' data(boston)
#' grb<-fgr1st(boston[,1:13])
#' grbu<-fundr(grb[[2]][,1:2])
fundr<-function(gr){
	gr<-matrix(gr,ncol=2)
	n<-length(gr[,1])
	for(i in 1:n) {
		if(gr[i,1]>gr[i,2]){jj<-gr[i,1]
			gr[i,1]<-gr[i,2]
			gr[i,2]<-jj
		}
	}
	gr0<-c(0,0)
	ind<-1:n
	rnk<-rank(gr[1:n,1],ties.method="first")
	ind[rnk]<-ind
	gr<-gr[ind,]
	jm<-max(gr[,1])
	for(j in 1:jm){
		ind1<-(1:n)[gr[,1]==j]
		li<-length(ind1)
		if(li>0){ind2<-ind1
			rnk<-rank(gr[ind1,2],ties.method="first")
			ind1[rnk]<-ind1
			gr[ind2,]<-gr[ind1,]
			gr2<-gr[ind2,2]
			gr2<-c(0,gr2)
			ad<-diff(gr2)
			lgr2<-length(ad)
			ind3<-(1:lgr2)[ad>0]
			gr1<-gr[ind2[ind3],]
			gr0<-rbind(gr0,gr1)
		}
	}
	li<-length(gr0[,1])
	gr<-gr0[2:li,]
	gr<-matrix(gr,ncol=2)
	list(gr)
}
	