#' Estimates the number of false positives by interpolating the results of simulations. These are valid for 50<= n <=5000, 25<= k <=50000, 1<= nu <= 10 and for p0=0.01 and p0=0.05. For values outside this range the numbers must be simulated.
#'
#' @param n  The dimension of dependent variable.
#' @param k  The number of covariates.
#' @param p0 Cut-off P-value.
#' @param nu Order statistic.
#' @param gr Logical, if TRUE then p0<-p0/k as is the default for graphs.
#' @param nufp Requires a data set of previous simulations nufp.rda
#' @param nsim Number of simulations. 
#' @param kmx Maximum number of selected covariates, must be larger than nu, for example nu+10
#' @param idum Seed for the random number generator.
#' @return enfp Estimated number of false positives
#' @return mnfp Mean number of false positives.
#' @return hist   Histogram of number of false positives.
#' @examples 
#' a<-fnfp(100,20,0.01,1:5,nufp,nsim=1000,kmx=10)
fnfp<-function(n,k,p0,nu,nufp,gr=F,nsim=0,kmx=0,idum=1){
	if(nsim>0){
		knu<-length(nu)
		if(gr){p0=p0/k}
		tmp<-.Fortran(
			"simords",
			as.integer(n),
			as.integer(k),
			double(n*k),
			double(n),
			double(n),
			as.double(nu),
			as.integer(nsim),
			double(knu*kmx),
			as.integer(idum),
			integer(k),
			as.integer(kmx),
			as.integer(knu),
			double(kmx),
			double(knu),
			as.double(p0)
			)
		rss<-tmp[[8]]
		dim(rss)<-c(knu,kmx)
		hist<-rss/nsim
		mn<-tmp[[14]]
		mn<-tmp[[14]]/nsim
		mn<-c(k,n,mn)
		enfp<--1
	}
	else{       
		mn<- -1
		hist<--1
		if((n<25)||(n>5000)||(k<25)||(k>50000))stop("(n,k) out of range")
		if((abs(p0-0.01)>0)&(abs(p0-0.05)>0))stop("p0 out of range")
		if(p0==0.01){
			if(!gr){efp<-nufp[1:70,]}
			else{efp<-nufp[71:140,]}
		}
		if(p0==0.05){
			if(!gr){efp<-nufp[141:210,]}
			else{efp<-nufp[211:280,]}
		}
		if(k>50000){
			if(n>=5000){fne<-efp[70,(1+nu)]}
			else{	
				if(n<=25){fne<-efp[71,3]}
				else{ j<-1
				while(j <6){
					if((n>=efp[(70+j),2])&(n<=efp[70+j+1,2])){jj<-j
						j<-6}
					j<-j+1
				}
				n1<-efp[(70+jj),2]
				n2<-efp[(70+jj+1),2]
				a<-(n-n1)/(n2-n1)
				b<-(n2-n)/(n2-n1)
				fne<-a*efp[(60+jj+1),(1+nu)]+b*efp[(60+jj),(1+nu)]
				}
			}
		}
		else{ 
			j<-1
			while(j <70){
				if((k<=efp[j+1,1])&(k>=efp[j,1])){jj<-j}
				j<-j+1
			}
			jj1<-jj-6
			jj2<-jj+1
			k1<-efp[jj1,1]
			k2<-efp[jj2,1]
#			print(c(jj,jj1,jj2,k1,k2))
			i<-1
			while(i<7){
#				print(c(n,efp[(jj1+i),2],efp[(jj1+i-1),2]))
				if((n<=efp[(jj1+i),2])&(n>=efp[(jj1+i-1),2])){ii<-i
					i<-6}
				i<-i+1
			}
			n1<-efp[(jj1+ii-1),2]
			n2<-efp[(jj1+ii),2]
			a<-(n-n1)/(n2-n1)
			b<-(n2-n)/(n2-n1)
			efp1<-a*efp[jj1+ii,]+b*efp[jj1+ii-1,]
			efp2<-a*efp[jj2+ii,]+b*efp[jj2+ii-1,]
			a<-(k-k1)/(k2-k1)
			b<-(k2-k)/(k2-k1)
			enfp<-a*efp2[2+nu]+b*efp1[2+nu]
		}
		if(gr){enfp<-k*enfp}
	
	}
		list(enfp,mn,hist)
}

