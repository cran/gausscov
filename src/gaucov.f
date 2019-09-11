C
C
C
C
      subroutine fststepwise(y,x,n,k,x2,res,ia,alpha,kmax,pp,kmax1,kexc
     $     ,intercept,nu,minss,chkintercept,misclass,mscls)
      integer n,k,kmax,kmax1
      double precision y(n),x(n,k),x2(n),res(n) ,pp(kmax1,2),minss(kmax1
     $     )
      integer ia(k+1),kexc(k+1),mscls(kmax1)
      logical intercept,chkintercept,misclass
      double precision alpha,nu
C
      integer icount,ks,ic,ik,kr,nex
      double precision ss0,ss1,amss1,pval,util1,util2,nx1,b,cf,pval1,mi
      double precision betai
C
C
      ic=0 
      do 1 j=1,k
         ia(j)=0
 1    continue
C
      if(intercept.and..not.chkintercept) then
         ia(1)=1
         b=0d0
         do 5 i=1,n
            b=b+y(i)
 5       continue
         b=b/dble(n)
         ss0=0d0
         do 6 i=1,n
            res(i)=y(i)-b
            y(i)=res(i)
            ss0=ss0+res(i)**2
 6       continue
         do 520 ik=2,k
            mi=0d0
            do 500 i=1,n
               mi=mi+x(i,ik)
 500        continue
            mi=mi/dble(n)
            do 510 i=1,n
               x(i,ik)=x(i,ik)-mi
 510        continue
 520     continue
         ks=1
         icount=1
      else
         ss0=0d0
         do 7 i=1,n
            ss0=ss0+y(i)**2
            res(i)=y(i)
 7       continue
         ks=0
         icount=0
      endif
C
      nex=0
      do 9 ik=1,k
         if(intercept.and.kexc(ik).ge.1) ia(kexc(ik)+1)=1
         if(.not.intercept.and.kexc(ik).ge.1) ia(kexc(ik))=1
         if(kexc(ik).gt.0) nex=nex+1
 9    continue
      kr=0
      do 19 ik=1,k
         if(ia(ik).eq.1) kr=kr+1
 19    continue
      kr=k-kr
C
 2    continue
      if(ks.eq.k) return
c
c
c
c
      ks=ks+1
      amss1=1d20
      do 40 kk=1,k
         if(ia(kk).eq.1) goto 40
         b=0d0
         nx1=0d0
         do 15 i=1,n
            b=b+x(i,kk)*res(i)
            nx1=nx1+x(i,kk)**2
 15      continue
         b=b/nx1
         ss1=0d0
         do 16 i=1,n
            ss1=ss1+(res(i)-b*x(i,kk))**2
 16      continue
         if(ss1.lt.amss1) then
            ic=kk
            amss1=ss1
            do 35 i=1,n
               x2(i)=x(i,kk)
 35         continue
         endif
 40   continue
      util1=1d0-amss1/ss0
c      if(intercept) then
c         util2=dble(n-icount-1)/2d0
c      else
         util2=dble(n-icount-1)/2d0   
c      endif   
      pval1=betai(util1,0.5d0,util2)
      pval=1d0-betai(pval1,dble(kr+2-icount)-nu,nu)
      if(pval.lt.alpha) then
         if(.not.intercept) icount=icount+1
         if(intercept.and.chkintercept) icount=icount+1
         minss(icount)=amss1
         b=0d0
         nx1=0d0
         do 45 i=1,n
            b=b+x2(i)*res(i)
            nx1=nx1+x2(i)**2
 45      continue
         b=b/nx1
         ss0=0d0
         cf=dsqrt(dble(n)/nx1)
         nx1=0d0
         do 50 i=1,n
            res(i)=res(i)-b*x2(i)
            ss0=ss0+res(i)**2
            x2(i)=cf*x2(i)
            nx1=nx1+x2(i)**2
 50      continue
         if(misclass) then
            mscls(icount) =0
            do 120 i=1,n
               mscls(icount)=mscls(icount)+min(1,iabs(idnint(res(i))))
 120        continue
         endif
         pp(icount,1)=dble(ic)
         if(intercept.and..not.chkintercept) pp(icount,1)=dble(ic)-1d0
         pp(icount,2)=pval
         if(icount.eq.kmax)  return
         if(intercept.and.icount+nex.eq.k-1) then
            kmax=icount
            return
         endif
         if(.not.intercept.and.icount+nex.eq.k) then
            kmax=icount
            return
         endif
         if(intercept.and..not.chkintercept) icount=icount+1
         ia(ic)=1
         do 60 kk=1,k
            if(ia(kk).eq.1) goto 60
            b=0d0
            do 61 i=1,n
c               b=b+x(i,kk)*xx(i,icount)
               b=b+x(i,kk)*x2(i)
 61         continue
            b=b/dble(n)
            do 62 i=1,n
               x(i,kk)=x(i,kk)-b*x2(i)
 62         continue
 60      continue
         goto 2
      else
         kmax=icount
         if(intercept.and..not.chkintercept) kmax=icount-1
         return
      endif
      end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       
C
C
      subroutine robstepwise(y,x,n,k,xx,x1,x2,beta,ia,alpha,kmax,pp
     $     ,beta0,cn,cpp,sig,res,res1,yy,kmax1,kexc,offset
     $     ,chkoffset,nu)
      integer n,k,kmax,kmax1
      double precision y(n),x(n,k),xx(n,kmax1),x1(n),x2(n) ,beta(kmax1)
     $     ,pp(kmax1,2) ,beta0(kmax1),res(n),res1(n),yy(n)
      integer ia(k),kexc(k)
      double precision alpha,cn,sig,cpp,nu
      logical offset,chkoffset
C
      integer icount,ks,ic,k1,nex
      double precision ss0,ss1,amss1,pval,sd10,sd20,sig1,b
     $     ,util1,sx,fc
      double precision ssrho(3)
      double precision rhoh,psih,psih1,betai
C
      ic=1 
      do 1 j=1,k
         ia(j)=0
 1    continue
      if(offset.and..not.chkoffset) then
         do 3 i=1,n
            xx(i,1)=1d0
 3       continue
         ia(1)=1
         k1=1
         call  orthrobreg(y,xx,yy,res,n,k1,beta,beta0,cn,sig,ssrho )
         ss0=ssrho(1)
         ks=1
         icount=1
      else
         ss0=0d0
         do 6 i=1,n
            ss0=ss0+rhoh(y(i)/sig,cn)
 6       continue
         ks=0
         icount=0
      endif
      nex=0
      do 9 ik=1,k
         if(offset.and.kexc(ik).ge.1) ia(kexc(ik)+1)
     $        =1
         if(.not.offset.and.kexc(ik).ge.1) ia(kexc(ik))=1
         if(kexc(ik).gt.0) nex=nex+1
 9    continue
      kr=0
      do 19 ik=1,k
         if(ia(ik).gt.0) kr=kr+1
 19    continue
      kr=k-kr
C
 2    continue
      if(ks.eq.k) return
      ks=ks+1
      amss1=1d20
C
      do 40 kk=1,k
         if(ia(kk).eq.1) goto 40
         do 11 i=1,n
            x1(i)=x(i,kk)
 11      continue
         do 14 j=1,icount
            b=0d0
            do 12 i=1,n
               b=b+x1(i)*xx(i,j)
 12         continue
            b=b/dble(n)
            do 13 i=1,n
               x1(i)=x1(i)-b*xx(i,j)
 13         continue
 14      continue
         sx=0d0
         do 15 i=1,n
            sx=sx+x1(i)**2
 15      continue
         if(sx.lt.1d-6) goto 40
         fc=dsqrt(dble(n)/sx)
         do 16 i=1,n
            xx(i,ks)=fc*x1(i)
 16      continue
         call orthrobreg(y,xx,yy,res,n,ks,beta,beta0,cn,sig,ssrho)
         ss1=ssrho(1)
         if(ss1.lt.amss1) then
            ic=kk
            amss1=ss1
            sd10=ssrho(2)
            sd20=ssrho(3)
            do 30 i=1,n
               res1(i)=res(i)
               x2(i)=fc*x1(i)
 30         continue
         endif
 40   continue
      util1=2d0*sd20*(ss0-amss1)/sd10
      pval=betai(util1/(4d2+util1),0.5d0,2d2)
      pval=1d0-betai(pval,dble(kr+2-icount)-nu,nu)
      if(pval.lt.alpha) then
         if(.not.offset) icount=icount+1
         if(offset.and.chkoffset) icount=icount+1
         pp(icount,1)=dble(ic)
         if(offset.and..not.chkoffset) pp(icount,1)=dble(ic)-1d0

         pp(icount,2)=pval
         if(icount.eq.kmax) return
         if(offset.and.icount+nex.eq.k-1) then
            kmax=icount
            return
         endif
         if(.not.offset.and.icount+nex.eq.k) then
            kmax=icount
            return
         endif
         if(offset.and..not.chkoffset) icount=icount+1
         ia(ic)=1
         do 45 i=1,n
            xx(i,icount)=x2(i)
 45      continue
         sig1=0d0
         do 50 i=1,n
            sig1=sig1+psih(res1(i)/sig,cn)**2
 50      continue
         sig1=dsqrt(sig1/(dble(n-ks)*cpp))*sig
         sig=sig1
         ss0=0d0
         sd10=0d0
         sd20=0d0
         do 60 i=1,n
            ss0=ss0+rhoh(res1(i)/sig,cn)
            sd10=sd10+psih(res1(i)/sig,cn)**2
            sd20=sd20+psih1(res1(i)/sig,cn)
 60      continue
         goto 2
      else
         kmax=icount
         if(offset.and..not.chkoffset) kmax=icount-1
         return
      endif
      end
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      subroutine robreg(y,x,yy,xx,xinv,n,k,d,r,beta,res,beta0,cn,sig
     $     ,ssrho,cbb)
      integer n,k
      double precision y(n),x(n,k),yy(n),xx(n,k),xinv(k,k),d(k),r(k)
     $     ,beta(k),res(n),beta0(k)
      double precision cn,sig,ssrho(3),cbb
      double precision sr1,sr2,ss,cn0,sig0
      double precision rhoh,psih,psih1
      integer ic
C
      logical inv
      inv=.false.
C
C
      do 10 i=1,n
         res(i)=y(i)
 10   continue
      call qrdecom(xx,n,k,d,r,inv)
      sig0=sig
      sr1=1d10
      ic=0
      cn0=cn
 100  continue
      ic=ic+1
      do 2 i=1,n
         yy(i)=psih(res(i)/sig0,cn0)*sig0
  2    continue
      call lsqqr(xx,yy,n,k,d,r,beta0,xinv,inv)
      do 35 j=1,k
         beta(j)=beta(j)+1.75d0*beta0(j)
 35   continue

      sr2=0d0
      do 50 i=1,n
         ss=0d0
         do 40 j=1,k
            ss=ss+x(i,j)*beta(j)
 40      continue
         res(i)=y(i)-ss
         sr2=sr2+rhoh(res(i)/sig0,cn0)
 50   continue
      if(sr1-sr2.gt.1d-2*sr2) then
         sr1=sr2
         goto 100
      endif
      sig=0d0
      do 55 i=1,n
         sig=sig+psih(res(i)/sig0,cn0)**2
 55   continue
      sig=dsqrt(sig/(cbb*dble(n-k)))*sig0
      if(dabs(sig/sig0-1d0).gt.1d-2) then
         sig0=sig
         goto100
      endif
      ssrho(1)=0d0
      ssrho(2)=0d0
      ssrho(3)=0d0
      do 60 i=1,n
         ssrho(1)=ssrho(1)+rhoh(res(i)/sig,cn)
         ssrho(2)=ssrho(2)+psih(res(i)/sig,cn)**2
         ssrho(3)=ssrho(3)+psih1(res(i)/sig,cn)
 60   continue
      return
      end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Calculates the least squares estimate for y=x*beta using the QR-
C     decomposition. The matrices x and y are overwritten. The coefficients
C     are returned in beta and their variances in d. Requires
C                              QR.f
C
      subroutine lsqqr(x,y,n,k,d,r,beta,x2inv,inv)
      integer n,k
      double precision x(n,k),y(n),d(k),r(k),beta(k),x2inv(k,k)
      logical inv
C
      double precision delta,sum
C
      do 30 j=1,k
         sum=0d0
         do 10 i=j,n
            sum=sum+x(i,j)*y(i)
 10      continue
         delta=sum/r(j)
         do 20 i=j,n
            y(i)=y(i)-delta*x(i,j)
 20      continue
 30   continue
C
      call qrsolv(x,y,n,k,d,beta)
      if(.not.inv) return
      do 60 j=1,k
         do 40 i=1,k
            y(i)=0d0
 40      continue
         y(j)=1d0
         call rsolv(x,y(1:k),n,k,d,r)
         do 50 i=1,k
           x2inv(i,j)=r(i)
 50      continue 
 60   continue
      do 80 j=1,k
         d(j)=0d0
         do 70 i=j,k
            d(j)=d(j)+x2inv(j,i)**2
 70      continue
 80   continue
      return
      end
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      subroutine qrdecom(x,m,n,d,r,sing)
      integer m,n
      double precision x(m,n),d(n),r(n)
      logical sing
C
      double precision beta,scl,sigma,sum
C
C
      sing=.false.
      do 80 j=1,n
         scl=0d0
         do 10 i=j,m
            scl=dmax1(scl,dabs(x(i,j)))
 10      continue
         if(scl.eq.0d0) then
            sing=.true.
            return
         endif
         do 20 i=j,m
            x(i,j)=x(i,j)/scl
 20      continue
         sum=0d0
         do 30 i=j,m
            sum=sum+x(i,j)**2
 30      continue
         sigma=dsign(dsqrt(sum),x(j,j))
         x(j,j)=x(j,j)+sigma
         r(j)=sigma*x(j,j)
         d(j)=-scl*sigma
         do 70 jj=j+1,n
            sum=0d0
            do 40 k=j,m
               sum=sum+x(k,jj)*x(k,j)
 40         continue
            beta=sum/r(j)
            do 60 k=j,m
               x(k,jj)=x(k,jj)-beta*x(k,j)
 60         continue
 70      continue
 80   continue
      return
      end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
C
      subroutine qrsolv(x,y,n,k,d,beta)
      integer n,k
      double precision x(n,k),y(n),d(k),beta(k)
C
      double precision sum
C
      beta(k)=y(k)/d(k)
C
      do 20 i=k-1,1,-1
         sum=0d0
         do 10 j=i+1,k
            sum=sum+x(i,j)*beta(j)
 10      continue
         beta(i)=(y(i)-sum)/d(i)
 20   continue
C
      return
      end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
      subroutine rsolv(x,y,n,k,d,beta)
      integer n,k
      double precision x(n,k),y(k),d(k),beta(k)
C
      double precision sum
C
      return
      beta(k)=y(k)/d(k)
C
      
      do 20 i=k-1,1,-1
         sum=0d0
         do 10 j=i+1,k
            sum=sum+x(i,j)*beta(j)
 10      continue
         beta(i)=(y(i)-sum)/d(i)
 20   continue
C
      return
      end
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      double precision function betai(x,a,b)
      double precision x,a,b
C
      double precision betacf, gammln
      double precision bt
C
C
C
      if(x.eq.0d0.or.x.eq.1d0) then
         bt=0d0
      else
         bt=dexp(gammln(a+b)-gammln(a)-gammln(b)+a*dlog(x)+b*dlog(1d0-x)
     $        )
      endif
      if(x.le.(a+1d0)/(a+b+2d0)) then
         betai=bt*betacf(a,b,x)/a
         return
      else
         betai=1d0-bt*betacf(b,a,1d0-x)/b
         return
      endif
C     
      end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      double precision function betacf(a,b,x)
      double precision a,b,x
C
      double precision eps,fpmin
      integer maxit
      parameter(maxit=100,eps=4d-10,fpmin=1d-20)
C
      double precision aa,c,d,del,h,qab,qam,qap
      integer m2,m
C
      qab=a+b
      qap=a+1d0
      qam=a-1d0
      c=1d0
      d=1d0-qab*x/qap
      if(dabs(d).lt.fpmin) d=fpmin
      d=1d0/d
      h=d
      do 11 m=1,maxit
         m2=2*m
         aa=dble(m)*(b-dble(m))*x/((qam+dble(m2))*(a+dble(m2)))
         d=1d0+aa*d
         if(dabs(d).lt.fpmin) d=fpmin
         c=1d0+aa/c
         if(dabs(c).lt.fpmin) c=fpmin
         d=1d0/d
         h=h*d*c
         aa=-(a+dble(m))*(qab+dble(m))*x/((a+dble(m2))*(qap+dble(m2)))
         d=1d0+aa*d
         if(dabs(d).lt.fpmin)d=fpmin
         c=1d0+aa/c
         if(dabs(c).lt.fpmin) c=fpmin
         d=1d0/d
         del=d*c
         h=h*del
         if(dabs(del-1d0).lt.eps) goto 1
 11   continue
C
 1    betacf=h
      return
      end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      double precision function gammln(xx)
      double precision xx
      integer j
      double precision ser,stp,tmp,x,y
      double precision cof(6)
      save cof,stp
      data cof,stp/76.18009172947146d0,-86.50532032941677d0,24
     $     .01409824083091d0,-1.231739572450155d0,0.1208650973866179d-2,
     $     -0.53952393849553d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5)*dlog(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
         y=y+1d0
         ser=ser+cof(j)/y
 11   continue
      gammln=tmp+dlog(stp*ser/x)
      return
      end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      subroutine edge(edg,ne,kmax,ke,ned)
      integer ne,kmax,ned
      integer edg(ne,2),ke(kmax)
C
      integer i0,i1,i,j
C
C
      call iquicksort(edg,ne,2,1)
C
      i0=1
      i1=1
      i2=1
      i=1
      j=edg(1,1)
      ic=0
 10   continue
      if(edg(i,1).eq.j) then
         ic=ic+1
         ke(ic)=edg(i,2)
         if(i.eq.ne) goto 15
         if(i.lt.ne) then
            i=i+1
            goto 10
         endif
      endif
 15   continue
      if(ic.gt.kmax) return
      call iquicksort(ke,ic,1,1)
      edg(i2,2)=ke(1)
      edg(i2,1)=j
      do 20 ii=2,ic
         if(ke(ii).eq.ke(ii-1)) goto 20
         i2=i2+1
         edg(i2,2)=ke(ii)
         edg(i2,1)=edg(i0,1)
 20   continue
      if(i.eq.ne) then
         ned=i2
         return
      endif
      i2=i2+1
      if(i.lt.ne) then
         i0=i
         j=edg(i,1)
         ic=0
         goto 10
      endif
      return
      end
C
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C
C
C
      subroutine graphst(xxx,x,n,k,y,x2,res,ia,alpha,kmax,pp,kmax1,grph
     $     ,ne,kexc,offset,nu,minss,chkoffset,mscls)
      integer n,k,kmax,kmax1,ne
      double precision xxx(n,k),x(n,k),y(n),x2(n),res(n),pp(kmax1,2)
     $     ,minss(kmax1)
      integer ia(k),grph(k*kmax1,2),kexc(k),mscls(kmax1)
      double precision alpha,nu
C
      logical offset,chkoffset,misclass
      integer ij,kmx
C
C
      misclass=.false.
      ne=0
      do 21 j=2,k
         do 6 j1=1,k
           do 5 i1=1,n
              x(i1,j1)=xxx(i1,j1)
 5         continue
 6      continue
         do 10 i=1,n
            y(i)=x(i,j)
 10      continue
         kexc(1)=j-1
         kmx=kmax
         call fststepwise(y,x,n,k,x2,res,ia,alpha,kmx,pp,kmax1,kexc
     $        ,offset,nu,minss,chkoffset,misclass,mscls)
         if(kmx.ge.1) then
            do 15, ij=1,kmx
               ne=ne+1
               if(j-1.lt.idnint(pp(ij,1))) then
                  grph(ne,1)=j-1
                  grph(ne,2)=idnint(pp(ij,1))
               else
                  grph(ne,1)=idnint(pp(ij,1))
                  grph(ne,2)=j-1
               endif
 15         continue
         endif
 21   continue
C
      return
      end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
C
      subroutine graphstst(xxx,x,n,k,y,x2,res,ia,alpha,kmax,pp,kmax1
     $     ,grph,ne,kexc,nedge,offset,nu,minss,chkoffset,mscls)
      integer n,k,kmax,kmax1,ne,nedge
      double precision x(n,k),y(n),x2(n),res(n),pp(kmax1,2),xxx(n,k)
     $     ,minss(kmax1)
      integer ia(k+1),grph(nedge,2),kexc(k+1),mscls(kmax1)
      double precision alpha,nu
      logical offset,chkoffset,misclass
C
      integer ij,kmx,ist
C
C
      misclass=.false.
      ne=0
      do 21,j=2,k

         do 10 i=1,n
            y(i)=xxx(i,j)
 10      continue
         do 11 iz=1,k
            kexc(iz)=0
 11      continue
         ist=1
         kexc(ist)=j-1
c
 12      continue
         do 6 j1=1,k
            do 5 i1=1,n
               x(i1,j1)=xxx(i1,j1)
 5          continue
 6       continue
         kmx=kmax
         call fststepwise(y,x,n,k,x2,res,ia,alpha,kmx,pp,kmax1,kexc
     $        ,offset,nu,minss,chkoffset,misclass,mscls)
      if(kmx.ge.1) then
         do 15, ij=1,kmx
            ne=ne+1
            ist=ist+1
            if(ne.gt.nedge) return
            kexc(ist)=idnint(pp(ij,1))
            if(j-1.lt.idnint(pp(ij,1))) then
               grph(ne,1)=j-1
               grph(ne,2)=idnint(pp(ij,1))
            else
               grph(ne,1)=idnint(pp(ij,1))
               grph(ne,2)=j-1
            endif
 15      continue
         goto 12
      endif
 21   continue
      return
      end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      SUBROUTINE QUICKSORT(A,N,SPALTE,INDEX)
C       Sortiert wird die Matrix A zeilenweise nach den Elementen der
C       Index-Spalte
      INTEGER N,T,NUK(17),NOK(17),SPALTE,INDEX
      DOUBLE PRECISION A(N,SPALTE)
      T=0
      NU=1
      NO=N
      K=0
1     IF (NU.LT.NO) THEN
            CALL TEILE(A,NU,NO,K,N,SPALTE,INDEX)
            T=T+1
            IF (K-NU.LT.NO-K) THEN
                NUK(T)=K+1
                NOK(T)=NO
                NO=K-1
            ELSE
                NUK(T)=NU
                NOK(T)=K-1
                NU=K+1
            ENDIF
            GOTO 1
      ENDIF
      IF (T.GT.0) THEN
            NU=NUK(T)
            NO=NOK(T)
            T=T-1
            GOTO 1
      ENDIF
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      SUBROUTINE TEILE(A,NU,NO,K,N,SPALTE,INDEX)
      INTEGER SPALTE
      DOUBLE PRECISION A(N,SPALTE)
      DOUBLE PRECISION H, HA(50)
      DO 5 K=1,SPALTE
5     HA(K)=A(NU,K)
      H=A(NU,INDEX)
      NUH=NU+1
      NOH=NO
10    DO 20 K=NOH,NUH,-1
      IF (A(K,INDEX).LT.H) THEN
        DO 15 J=1,SPALTE
15      A(NUH-1,J)=A(K,J)
        NOH=K-1
        GOTO 30
      ENDIF
20    CONTINUE
      DO 25 I=1,SPALTE
25    A(K,I)=HA(I)
      RETURN
30    DO 40 K=NUH,NOH
      IF (A(K,INDEX).GT.H) THEN
        DO 35 J=1,SPALTE
35      A(NOH+1,J)=A(K,J)
        NUH=K+1
        GOTO 10
      ENDIF
40    CONTINUE
      DO 45 I=1,SPALTE
45    A(K,I)=HA(I)
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      SUBROUTINE IQUICKSORT(A,N,SPALTE,INDEX)
C       Sortiert wird die Matrix A zeilenweise nach den Elementen der
C       Index-Spalte
      INTEGER N,T,NUK(17),NOK(17),SPALTE,INDEX
      INTEGER  A(N,SPALTE)
      T=0
      NU=1
      NO=N
      K=0
1     IF (NU.LT.NO) THEN
            CALL ITEILE(A,NU,NO,K,N,SPALTE,INDEX)
            T=T+1
            IF (K-NU.LT.NO-K) THEN
                NUK(T)=K+1
                NOK(T)=NO
                NO=K-1
            ELSE
                NUK(T)=NU
                NOK(T)=K-1
                NU=K+1
            ENDIF
            GOTO 1
      ENDIF
      IF (T.GT.0) THEN
            NU=NUK(T)
            NO=NOK(T)
            T=T-1
            GOTO 1
      ENDIF
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      SUBROUTINE ITEILE(A,NU,NO,K,N,SPALTE,INDEX)
      INTEGER SPALTE
      INTEGER A(N,SPALTE)
      INTEGER H, HA(50)
C      
C
      DO 5 K=1,SPALTE
5     HA(K)=A(NU,K)
      H=A(NU,INDEX)
      NUH=NU+1
      NOH=NO
10    DO 20 K=NOH,NUH,-1
      IF (A(K,INDEX).LT.H) THEN
        DO 15 J=1,SPALTE
15      A(NUH-1,J)=A(K,J)
        NOH=K-1
        GOTO 30
      ENDIF
20    CONTINUE
      DO 25 I=1,SPALTE
25    A(K,I)=HA(I)
      RETURN
30    DO 40 K=NUH,NOH
      IF (A(K,INDEX).GT.H) THEN
        DO 35 J=1,SPALTE
35      A(NOH+1,J)=A(K,J)
        NUH=K+1
        GOTO 10
      ENDIF
40    CONTINUE
      DO 45 I=1,SPALTE
45    A(K,I)=HA(I)
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      double precision function psih1(x,cn)
      double precision x,cn
C
c      double precision u
C
c      if(dabs(cn*x).ge.1.5d1) then
c         psih1=0d0
c      else
c         u=dexp(cn*x)
c         psih1=2d0*cn*u/(1d0+u)**2
c      endif
      if(dabs(x).le.cn) then
         psih1=1d0
      else
         psih1=0d0
      endif
      return
      end
C
C
C
C
      double precision function psih(x,cn)
      double precision x,cn
C
      double precision dsign
c      double precision u
C
c      if(dabs(cn*x).le.1.5d1) then
c         u=dexp(cn*x)
c         psih=(u-1d0)/(u+1d0)
c      else
c        psih=dsign(1d0,x)
c      endif
      if(dabs(x).le.cn) then
         psih=x
      else
         psih=cn*dsign(1d0,x)
      endif
      return
      end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
      subroutine orthrobreg(y,x,yy,res,n,ks,beta,beta0,cn,sig,ssrho)
      integer n,ks
      double precision y(n),x(n,ks),yy(n),res(n),beta(ks),beta0(ks)
      double precision cn,sig,ssrho(3)
      double precision sr1,sr2,ss,cbb,syx
      double precision rhoh,psih,psih1
      integer ic
C
C
      cbb=0.5836578d0
C
      do 1 i=1,n
         res(i)=y(i)
 1    continue
      sr1=1d10
      ic=0
 100  continue
      ic=ic+1
      do 2 i=1,n
         yy(i)=psih(res(i)/sig,cn)*sig
  2    continue
       do 20 j=1,ks
          syx=0d0
          do 10 i=1,n
             syx=syx+yy(i)*x(i,j)
 10       continue
          beta0(j)=syx/dble(n)
 20    continue
      do 35 j=1,ks
         beta(j)=beta(j)+1.75d0*beta0(j)
 35   continue
      sr2=0d0
      do 50 i=1,n
         ss=0d0
         do 40 j=1,ks
            ss=ss+x(i,j)*beta(j)
 40      continue
         res(i)=y(i)-ss
         sr2=sr2+rhoh(res(i)/sig,cn)
 50   continue
      if(sr1-sr2.gt.1d-3*sr2) then
         sr1=sr2
         goto 100
      endif
      ssrho(1)=0d0
      ssrho(2)=0d0
      ssrho(3)=0d0
      do 60 i=1,n
         ssrho(1)=ssrho(1)+rhoh(res(i)/sig,cn)
         ssrho(2)=ssrho(2)+psih(res(i)/sig,cn)**2
         ssrho(3)=ssrho(3)+psih1(res(i)/sig,cn)
 60   continue
      return
      end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      double precision function rhoh(x,cn)
      double precision x,cn
C
C
c      if(dabs(cn*x).ge.15d1) then
c         rhoh=dabs(x)
c      else
c        rhoh=2d0*dlog(0.5d0+0.5d0*dexp(cn*x))/cn-x
c      endif
      if(dabs(x).le.cn) then
         rhoh=0.5d0*x**2
      else
         rhoh=cn*dabs(x)-0.5d0*cn**2
      endif
      return
      end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
C
      subroutine lmmdch(y,x,n,k,xx,xxx,y1,y2,d,r,beta,xinv,ia,intercept
     $     ,ss,nv,ssr,p0)
      integer n,k
      double precision y(n),x(n,k),xx(n*k),xxx(n*k),y1(n),y2(n) ,d(k)
     $     ,r(k),beta(k),xinv(k**2),ss(2**k),ssr(2**k)
      double precision p0
      integer ia(k),nv(2**k,2)
      logical intercept
C
      double precision ss0,ss1,pval,util1,util2
      double precision betai
      integer id,ic,ks,ns,np,ni
      logical inv
C
      inv=.false.
      id=1
      ss(1)=0d0
      do 300 i=1,n
         ss(1)=ss(1)+y(i)**2
 300  continue
C
C
      ni=0
      do 40 iv=1,2**k-1
         ssr(iv)=0d0
         call decode(iv,k,ia)
         ks=0
         do 20 i=1,k
            ks=ks+ia(i)
 20     continue
        call xsubset1(x,xx,n,k,ks,ia,id)
        call lsq(xx,y,xxx,y1,n,ks,d,r,beta,xinv,y2,inv)
         ss(iv+1)=0d0
         do 30 i=1,n
            ss(iv+1)=ss(iv+1)+y2(i)**2
 30      continue
 40   continue
C
      ic=0
      do 80 iv =1,2**k-1
         call decode(iv,k,ia)
         ks=0
         do 45 i=1,k
            ks=ks+ia(i)
 45      continue
         np=0
         do 70 is=1,k
             if(ia(is).eq.1) then
               ia(is)=0
               ns=1
               do 50 iu=1,k
                  ns=ns+ia(iu)*2**(iu-1)
 50            continue
               ia(is)=1
               ss1=ss(iv+1)
               ss0=ss(ns)

               util1=1d0-ss1/ss0
               if(intercept) then
                  util2=dble(n-ks-1)/2d0
               else
                  util2=dble(n-ks)/2d0   
               endif   
               pval=1d0-betai(util1,0.5d0,util2)
               if(pval.gt.p0) goto 80
               np=np+1
            endif
 70      continue
           ni=ni+1
           nv(ni,1)=iv
           nv(ni,2)=np
           ssr(ni)=ss(iv+1)
 80      continue
         return
C
      end
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      subroutine decode(j,k,set)
      integer j,k
      integer set(k),jj

      jj=j
      do 10 i=1,k
         set(i)=0
 10   continue
      if(j.eq.0) return
      do 20 i =k-1,0,-1
         if(jj.ge.2**i) then
            set(i+1)=1
            jj=jj-2**i
         endif
 20   continue
      return
      end
C
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
      subroutine lsq(x,y,xx,yy,n,k,d,r,beta,x2inv,res,inv)
      integer n,k
      double precision x(n,k),y(n),xx(n,k),yy(n),d(k),r(k),beta(k)
     $     ,x2inv(k,k),res(n)
      logical inv
C
      double precision b
C
      do 20 i=1,n
         yy(i)=y(i)
         do 10 j=1,k
            xx(i,j)=x(i,j)
 10      continue
 20   continue
C
      call qrdecom(xx,n,k,d,r,inv)
      if(inv) return
      call  lsqqr(xx,yy,n,k,d,r,beta,x2inv,inv)
C
      do 40 i=1,n
         b=0d0
         do 30 j=1,k
            b=b+x(i,j)*beta(j)
 30      continue
         res(i)=y(i)-b
 40   continue
C
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      subroutine xsubset1(x,xx,n,k,ks,ia,id)
      integer n,k,ks,id
      double precision x(n,k),xx(n,ks)
      integer ia(k)
C
      integer ic
C
      ic=0
      do 20 j=1,k
         if(ia(j).eq.id) then
            ic=ic+1
            do 10 i=1,n
               xx(i,ic)=x(i,j)
 10         continue
         endif
 20   continue
      end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      subroutine genint(x,xx,n,k,kk,ord,ind)
      integer n,k,kk,ord
      integer ind(ord)
      double precision x(n,k),xx(n,kk)
C
C
      do 10 i=1,ord
         ind(i)=1
 10   continue
C
      ic=0
      do 40 j=1,kk
         do 30 i=1,n
            xx(i,j)=1d0
            do 20 jj=1,ord
               xx(i,j)=xx(i,j)*x(i,ind(jj))
 20         continue
 30      continue
         ic=ic+1
         call inact(ind,k,ord)
 40   continue
      return
      end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      subroutine degenint(ic,m,k,ord,ind,iut)
      integer m,k,ord
      integer ind(m,ord),ic(m),iut(ord)
C
C
      do 20 j=1,m
         call degenint1(ic(j),k,ord,iut)
         do 10 i=1,ord
            ind(j,i)=iut(i)-1
 10      continue
 20   continue
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      subroutine degenint1(ic0,k,ord,ind)
      integer ic0,k,ord
      integer ind(ord)
C
C
      do 10 i=1,ord
         ind(i)=1
 10   continue
C
      ic=0
 40   continue
      ic=ic+1
      if(ic.eq.ic0) return
      call inact(ind,k,ord)
      goto 40
      return
      end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      subroutine inact(ina,k,ord)
      integer ord,k
      integer ina(ord)
C
      integer j
c
      if(ina(ord).lt.k) then
         ina(ord)=ina(ord)+1
         return
      endif
      j=ord-1
 10   continue
      if(ina(j).lt.k) then
         ina(j)=ina(j)+1
         do 20 jj=j+1,ord
            ina(jj)=ina(j)
 20      continue
         return
      elseif(j.gt.1) then
         j=j-1
         goto 10
      endif
      return
      end
C
