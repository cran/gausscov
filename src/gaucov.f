C
C
C
      subroutine fstepwise(y,x,n,k,x2,res,ia,alpha,kmax,pp,kmax1,kexc
     $     ,intercept,nu,minss,ss01)
      integer n,k,kmax,kmax1
      double precision y(n),x(n,k),x2(n),res(n) ,pp(kmax1,2),minss(kmax1
     $     ),ss01(k+1)
      integer ia(k+1),kexc(k+1)
      logical intercept
      double precision alpha,nu
C
      integer icount,ks,ic,ik,kr,nex
      double precision ss0,ss1,amss1,pval,util1,util2,nx1,b,cf,pval1,mi
     $     ,ssy
      double precision betai
C
      ic=0 
      do 1 j=1,k+1
         ia(j)=0
 1    continue
C
      ssy=0d0
      do 3 i=1,n
         ssy=ssy+y(i)**2
 3    continue
       if(intercept) then
         ia(k)=1
         b=0d0
         do 5 i=1,n
            b=b+y(i)
 5       continue
         b=b/dble(n)
         ss0=0d0
         ss1=0d0
         do 6 i=1,n
           ss0=ss0+y(i)**2
            res(i)=y(i)-b
            ss1=ss1+res(i)**2
 6       continue
         util1=ss1/ss0
         pval=betai(util1,0.5d0*dble(n-1),0.5d0)
         pp(1,1)=dble(k)
         pp(1,2)=pval
         minss(1)=ss1
         ss01(1)=ss1/ss0
         do 520 ik=1,k-1
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
         ss0=ss1
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
c         if(intercept.and.kexc(ik).ge.1)then
c            kexc(ik)=kexc(ik)+1
c         endif
         if(kexc(ik).gt.0) then
            ia(kexc(ik))=1
            nex=nex+1
        endif
 9    continue
      kr=0
      do 19 ik=1,k
         if(ia(ik).eq.1) kr=kr+1
 19    continue
      kr=k-kr
C
 2    continue
      if(ks.eq.k) goto 600
C
C
C
      ks=ks+1
      amss1=ss0
      do 40 kk=1,k
         if(ia(kk).eq.1) goto 40
         b=0d0
         nx1=0d0
         do 15 i=1,n
            b=b+x(i,kk)*res(i)
            nx1=nx1+x(i,kk)**2
 15      continue
         if(nx1.lt.1d-6) goto 40
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
      ks=icount+1
      util1=amss1/ss0
      util2=dble(n-icount-1)/2d0   
      pval1=betai(util1,util2,0.5d0)
      pval=betai(pval1,nu,dble(kr+2-icount)-nu)
c      write(*,*) pval,alpha
      if(pval.lt.alpha) then
          icount=icount+1
         minss(icount)=amss1
         ss01(icount)=amss1/ss0
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
         pp(icount,1)=dble(ic)
c         if(intercept) pp(icount,1)=dble(ic)-1d0
         pp(icount,2)=pval
         if(icount.eq.kmax)  goto 600
         if(icount+nex.eq.k) then
            kmax=icount
            goto 600
         endif
         ia(ic)=1
         do 60 kk=1,k
            if(ia(kk).eq.1) goto 60
            b=0d0
            do 61 i=1,n
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
         goto 600
      endif
 600  continue
      return
      end
C
C
C      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       
C
C
      subroutine robstepwise(y,x,n,k,xx,x1,x2,beta,ia,alpha,kmax,pp
     $     ,beta0,cn,cpp,sig,res,res1,yy,kmax1,kexc,offset,nu,ssg,red,
     $     cnr)
      integer n,k,kmax,kmax1
      double precision y(n),x(n,k),xx(n,kmax1),x1(n),x2(n) ,beta(kmax1)
     $     ,pp(kmax1,2) ,beta0(kmax1),res(n),res1(n),yy(n),ssg(kmax1),
     $     cnr(3)
      integer ia(k+1),kexc(k+1)
      double precision alpha,cn,sig,cpp,nu
      logical offset,red
C
      integer icount,ks,ic,nex,id
      double precision ss0,ss1,amss1,pval,sd10,sd20,b
     $     ,util1,sx,fc
      double precision ssrho(3)
      double precision betai,sig0,rhoh,rrhoh
      logical scale
C
C
      id=1
      sd10=1d0
      sd20=1d0
      ic=1 
      do 1 j=1,k
         ia(j)=0
 1    continue
      if(offset) then
         ss0=0d0
         do 3 i=1,n
            xx(i,1)=1d0
            if(red) then
               ss0=ss0+rrhoh(y(i)/sig,cnr)
            else
               ss0=ss0+rhoh(y(i)/sig,cn)
            endif
 3       continue
         ia(k)=1
         ks=1
         scale=.true.
         call  orthrobreg(y,xx,yy,res,n,ks,beta,beta0,cn,cpp,sig,ssrho
     $        ,scale,red,cnr)
         ss1=ssrho(1)
         sd10=ssrho(2)
         sd20=ssrho(3)
         util1=2d0*sd20*(ss0-ss1)/sd10
         pval=betai(util1/(4d2+util1),0.5d0,2d2)
         pp(1,1)=dble(k)
         pp(1,2)=1d0-pval
         icount=1
         ssg(1)=sig
         ks=1
         ss0=ss1
         icount=1
      else
         ss0=0d0
         do 6 i=1,n
            if(red) then
               ss0=ss0+rrhoh(y(i)/sig,cnr)
            else
               ss0=ss0+rhoh(y(i)/sig,cn)
            endif
 6       continue
         ks=0
         icount=0
      endif
      nex=0
      do 9 ik=1,k
         if(kexc(ik).gt.0) then
            ia(kexc(ik))=1
            nex=nex+1
         endif
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
c         sig=ssg(kk+1)
         scale=.false.
         call orthrobreg(y,xx,yy,res,n,ks,beta,beta0,cn,cpp,sig,ssrho
     $        ,scale,red,cnr)
         ss1=ssrho(1)
         if(ss1.lt.amss1) then
            ic=kk
            amss1=ss1
            sd10=ssrho(2)
            sd20=ssrho(3)
            sig0=sig
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
         icount=icount+1
         pp(icount,1)=dble(ic)
c         if(offset) pp(icount,1)=dble(ic)-1d0
         pp(icount,2)=pval
         ia(ic)=1
         do 45 i=1,n
            xx(i,icount)=x2(i)
 45      continue
         ks=icount
         scale=.true.
         call orthrobreg(y,xx,yy,res,n,ks,beta,beta0,cn,cpp,sig,ssrho
     $        ,scale,red,cnr)
        ssg(icount)=sig
        if(icount.eq.kmax) return
        if(icount+nex.eq.k) then
           kmax=icount
           return
        endif
        ss0=ssrho(1)
        sd10=ssrho(2)
        sd20=ssrho(3)
        goto 2
      else
         kmax=icount
         return
      endif
      end
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
      subroutine  robregp(y,x,yy,xx,xxx,xinv,n,k,d,r,beta,res,beta0,cn
     $     ,sig,ssrho,cpp,ia,pp,q,intercept,res1,scale,red,cnr)
      integer n,k,q
      double precision y(n),x(n,k),yy(n),xx(n,k),xxx(n*k),xinv(k**2),d(k
     $     ),r(k),beta(k),res(n),beta0(k),pp(k,2),res1(n)
      double precision cn,sig,ssrho(3),cpp,cnr(3)
      integer ia(k)
      logical intercept,scale,red
C
      integer id
      double precision ss0,ss11,ss1,sd10,sd20,pval,pval1,util1,sig0
      double precision betat(1000),ssrho1(3)
      double precision betai,rhoh,psih,psih1,rrhoh,rpsih,rpsih1
C
      sig0=sig
      kk=k
      id=1
      do 20 i=1,k
         ia(i)=1
 20   continue
      ks=k
      call xsubset1(x,xx,n,k,ks,ia,id)
      scale=.true.
      sig=sig0  
      call robreg(y,xx,yy,xxx,xinv,n,ks,d,r,betat,res1,beta0,cn,sig ,
     $     ssrho,cpp,scale,red,cnr)
      ss1=ssrho(1)
      if(k.eq.1) then
         ss11=0d0
         ssrho1(1)=0d0
         ssrho1(2)=0d0
         ssrho1(3)=0d0
         do 100 i=1,n
            if(red) then
               ss11=ss11+rrhoh(res1(i)/sig0,cnr)
               ssrho1(1)=ssrho1(1)+rrhoh(y(i)/sig0,cnr)
               ssrho1(2)=ssrho1(2)+rpsih(y(i)/sig0,cnr)**2
               ssrho1(3)=ssrho1(3)+rpsih1(y(i)/sig0,cnr)
            else
               ss11=ss11+rhoh(res1(i)/sig0,cn)
               ssrho1(1)=ssrho1(1)+rhoh(y(i)/sig0,cn)
               ssrho1(2)=ssrho1(2)+psih(y(i)/sig0,cn)**2
               ssrho1(3)=ssrho1(3)+psih1(y(i)/sig0,cn)
            endif
 100     continue
         sd10=ssrho1(2)
         sd20=ssrho1(3)
         ss0=ssrho1(1)
         if(ss0.lt.ss11+1d-6) then
            pval=1d0
            pval1=1d0
         else
            util1=2d0*sd20*(ss0-ss11)/sd10
            pval1=1d0-betai(util1/(4d2+util1),0.5d0,2d2) 
         endif
         if(intercept.and.kk.eq.k) then
             pval=pval1
         else 
            pval=betai(pval1,1d0,dble(q+1-k))
         endif
          pp(1,1)=pval
          pp(1,2)=pval
          beta(1)=betat(1)
          return
       endif
       ks=k-1
       scale=.false.
c
       do 30 ik=1,k
          ia(ik)=0
          call xsubset1(x,xx,n,k,ks,ia,id)
          ia(ik)=1
          call robreg(y,xx,yy,xxx,xinv,n,ks,d,r,beta,res,beta0,cn,sig
     $        ,ssrho1,cpp,scale,red,cnr)
          sd10=ssrho1(2)
          sd20=ssrho1(3)
          ss0=ssrho1(1)

          if(ss0.lt.ss1+1d-6) then
             pval=1d0
             pval1=1d0
          else
             util1=2d0*sd20*(ss0-ss1)/sd10
             pval1=1d0-betai(util1/(4d2+util1),0.5d0,2d2)
             if(intercept.and.ik.eq.k) then
                pval=pval1
             else
              pval=betai(pval1,1d0,dble(q+1-k))
            endif
          endif
          pp(ik,1)=pval
          pp(ik,2)=pval1
   30   continue
        do 400 i=1,k
           beta(i)=betat(i)
 400    continue
      end
C
C

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      subroutine robreg(y,x,yy,xx,xinv,n,k,d,r,beta,res,beta0,cn,sig
     $     ,ssrho,cpp,scale,red,cnr)
      integer n,k
      double precision y(n),x(n,k),yy(n),xx(n,k),xinv(k,k),d(k),r(k)
     $     ,beta(k),res(n),beta0(k)
      double precision cn,sig,ssrho(3),cpp,cnr(3)
      double precision sr1,sr2,ss,sig0
      double precision rhoh,psih,psih1,rrhoh,rpsih,rpsih1
      integer ic
      logical scale,red
C
      logical inv
      inv=.false.
C
      do 20 i=1,n
         yy(i)=y(i)
         do 15 j=1,k
            xx(i,j)=x(i,j)
 15      continue
 20   continue
C
      do 1 i=1,k
         beta(i)=0d0
 1    continue
      do 10 i=1,n
         res(i)=y(i)
 10   continue
      call qrdecom(xx,n,k,d,r,inv)
      sig0=sig
      sr1=1d10
      ic=0
 100  continue
      ic=ic+1
      do 2 i=1,n
         if(red) then
            yy(i)=rpsih(res(i)/sig0,cnr)*sig0
         else
            yy(i)=psih(res(i)/sig0,cn)*sig0
         endif
 2    continue
       call lsqqr(xx,yy,n,k,d,r,beta0,xinv,inv)
 
      do 35 j=1,k
         beta(j)=beta(j)+1.5d0*beta0(j)
 35   continue
        sr2=0d0
      do 50 i=1,n
         ss=0d0
         do 40 j=1,k
            ss=ss+x(i,j)*beta(j)
 40      continue
         res(i)=y(i)-ss
         if(red) then
            sr2=sr2+rrhoh(res(i)/sig0,cnr)
         else
            sr2=sr2+rhoh(res(i)/sig0,cn)
         endif
 50   continue
      if(sr1-sr2.gt.1d-3*sr2) then
         sr1=sr2
         goto 100
      endif
      if(.not.scale) goto 200
      sig=0d0
      do 55 i=1,n
         if(red) then
            sig=sig+rpsih(res(i)/sig0,cnr)**2
         else
            sig=sig+psih(res(i)/sig0,cn)**2
         endif
 55   continue
       sig=dsqrt(sig/(cpp*dble(n-k)))*sig0
       if(dabs(sig/sig0-1d0).gt.1d-2) then
         sig0=sig
         goto 100
      endif
 200  continue
      ssrho(1)=0d0
      ssrho(2)=0d0
      ssrho(3)=0d0
      do 60 i=1,n
         if(red) then
            ssrho(1)=ssrho(1)+rrhoh(res(i)/sig,cnr)
            ssrho(2)=ssrho(2)+rpsih(res(i)/sig,cnr)**2
            ssrho(3)=ssrho(3)+rpsih1(res(i)/sig,cnr)
         else
            ssrho(1)=ssrho(1)+rhoh(res(i)/sig,cn)
            ssrho(2)=ssrho(2)+psih(res(i)/sig,cn)**2
            ssrho(3)=ssrho(3)+psih1(res(i)/sig,cn)
         endif
 60   continue
      return
      end
C
C
C
C
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
      parameter(maxit=200,eps=4d-20,fpmin=1d-20)
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
      integer i0,i1,i,j,ll
C
C
      do 1 i=1,ne
         if(edg(i,1).gt.edg(i,2)) then
            ll=edg(i,1)
            edg(i,1)=edg(i,2)
            edg(i,2)=ll
         endif
 1      continue
         
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
     $     ,ne,kexc,offset,nu,minss,nedge,ss01,kk)
      integer n,k,kmax,kmax1,ne,nedge,kk
      double precision xxx(n,k),x(n,kk),y(n),x2(n),res(n),pp(kmax1,2)
     $     ,minss(kmax1),ss01(kk)
      integer ia(kk+1),grph(k*kmax1,2),kexc(kk+1)
      double precision alpha,nu
C
      logical offset
      integer ij,kmx,istrt
C
C
      ne=0
      do 6 j=1,k
         do 4 j1=1,k
            do 3 ii=1,n
                  x(ii,j1)=xxx(ii,j1)
 3          continue
 4       continue
         do 10 i=1,n
            y(i)=xxx(i,j)
 10      continue
         kmx=kmax
         kexc(1)=j
         call fstepwise(y,x,n,kk,x2,res,ia,alpha,kmx,pp,kmax1,kexc
     $        ,offset,nu,minss,ss01)
         if(kmx.eq.0) goto 6
         if(kmx.eq.1.and.idnint(pp(1,1)).eq.0) goto 6
         if(idnint(pp(1,1)).eq.0) then
            istrt=2
         else
            istrt=1
         endif
         do 15, ij=2,kmx
            if(idnint(pp(ij,1)).ge.1) then
               ne=ne+1
               grph(ne,1)=j
               grph(ne,2)=idnint(pp(ij,1))
               if(ne.ge.nedge) return
            endif
 15      continue
 6    continue
C
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      subroutine graphstst(xxx,x,n,k,y,x2,res,ia,alpha,kmax,pp,kmax1
     $     ,grph,ne,kexc,nedge,offset,nu,minss,ss01,kk)
      integer n,k,kmax,kmax1,ne,nedge
      double precision x(n,kk),y(n),x2(n),res(n),pp(kmax1,2),xxx(n,k)
     $     ,minss(kmax1),ss01(kk)
      integer ia(kk+1),grph(k*kmax1,2),kexc(kk+1)
      double precision alpha,nu
      logical offset
C
      integer ij,kmx,ek,la
C
C
      ne=0
      do 100 j=1,k
         ek=1
         la=0
         do 1 j1=1,k
            kexc(j1)=0
 1       continue
         do 2 i=1,n
            y(i)=xxx(i,j)
 2       continue
 3       continue
         kexc(1)=j
         do 5 j1=1,k
            do 4 ii=1,n
c               if(offset) then
c                  x(ii,j1+1)=xxx(ii,j1)
c               else
                  x(ii,j1)=xxx(ii,j1)
c               endif
 4          continue
 5       continue
         kmx=kmax-ek
         call fstepwise(y,x,n,kk,x2,res,ia,alpha,kmx,pp,kmax1,kexc
     $        ,offset,nu,minss,ss01)
         if(kmx.le.1) goto 100
         la=la+1
         do 15, ij=2,kmx
            if(idnint(pp(ij,1)).ge.1) then
               ne=ne+1
               ek=ek+1
               kexc(ek)=idnint(pp(ij,1))
               grph(ne,1)=j
               grph(ne,2)=idnint(pp(ij,1))
               if(offset) then
                  if(grph(ne,2).eq.kk) grph(ne,2)=0
               endif
               if(ne.ge.nedge) return
            endif
 15      continue
         if(ek.lt.kmax) goto 3
 100  continue
C
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
         HA(K)=A(NU,K)
 5    CONTINUE
      H=A(NU,INDEX)
      NUH=NU+1
      NOH=NO
 10   CONTINUE
      DO 20 K=NOH,NUH,-1
         IF (A(K,INDEX).LT.H) THEN
            DO 15 J=1,SPALTE
               A(NUH-1,J)=A(K,J)
 15         CONTINUE
            NOH=K-1
            GOTO 30
         ENDIF
20    CONTINUE
      DO 25 I=1,SPALTE
         A(K,I)=HA(I)
 25   CONTINUE
      RETURN
 30   CONTINUE
      DO 40 K=NUH,NOH
         IF (A(K,INDEX).GT.H) THEN
            DO 35 J=1,SPALTE
               A(NOH+1,J)=A(K,J)
 35         CONTINUE
            NUH=K+1
            GOTO 10
         ENDIF
40    CONTINUE
      DO 45 I=1,SPALTE
         A(K,I)=HA(I)
 45   CONTINUE
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
         HA(K)=A(NU,K)
 5    CONTINUE
      H=A(NU,INDEX)
      NUH=NU+1
      NOH=NO
10    DO 20 K=NOH,NUH,-1
      IF (A(K,INDEX).LT.H) THEN
        DO 15 J=1,SPALTE
           A(NUH-1,J)=A(K,J)
 15     CONTINUE
        NOH=K-1
        GOTO 30
      ENDIF
20    CONTINUE
      DO 25 I=1,SPALTE
         A(K,I)=HA(I)
 25   CONTINUE
      RETURN
30    DO 40 K=NUH,NOH
      IF (A(K,INDEX).GT.H) THEN
        DO 35 J=1,SPALTE
           A(NOH+1,J)=A(K,J)
 35     CONTINUE
        NUH=K+1
        GOTO 10
      ENDIF
40    CONTINUE
      DO 45 I=1,SPALTE
         A(K,I)=HA(I)
 45   CONTINUE
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
c      return
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      double precision function rpsih1(x,cn)
      double precision x,cn(3),ax
C
      double precision dsign
C
      ax=dabs(x)
      if(ax.le.cn(1)) then
         rpsih1=1d0
      elseif(ax.ge.cn(1).and.ax.le.cn(2)) then
         rpsih1=0d0      
      elseif(cn(2).le.ax.and.ax.le.cn(3)) then
         rpsih1=dsign(1d0,x)*cn(1)/(cn(3)-cn(2))              
      else
         rpsih1=0d0
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
c      return
      if(dabs(x).le.cn) then
         psih=x
      else
         psih=cn*dsign(1d0,x)
      endif
      return
      end
C
C
C
      double precision function rpsih(x,cn)
      double precision x,cn(3),ax
C
      double precision dsign
c      double precision u
C
      ax=dabs(x)
      if(ax.le.cn(1)) then
         rpsih=x
      elseif(ax.ge.cn(1).and.ax.le .cn(2)) then
         rpsih=cn(1)*dsign(1d0,x)
      elseif(ax.ge.cn(2).and.ax.le.cn(3)) Then
         rpsih=dsign(1d0,x)*cn(1)*(cn(3)-ax)/(cn(3)-cn(2))
      else
         rpsih=0d0
      endif
      return
      end
C
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
c      return
      if(dabs(x).le.cn) then
         rhoh=0.5d0*x**2
      else
         rhoh=cn*dabs(x)-0.5d0*cn**2
      endif
      return
      end
C
C
C
C
      double precision function rrhoh(x,cn)
      double precision x,cn(3),ax
C
C
      rrhoh=0d0
      ax=dabs(x)
      if(ax.le.cn(1)) then
         rrhoh=0.5d0*ax**2
      elseif(ax.ge.cn(1).and.ax.le.cn(2)) then
         rrhoh=cn(1)*ax-0.5d0*cn(1)**2

      elseif(ax.ge.cn(2).and.ax.le.cn(3)) then
         rrhoh=cn(1)*cn(2)-0.5d0*cn(1)**2+0.5d0*cn(1)*((cn(3)-cn(2))**2-
     $  (cn(3)-ax)**2)/(cn(3)-cn(2))
      elseif(ax.ge.cn(3)) then
         rrhoh=cn(1)*cn(2)-0.5d0*cn(1)**2+0.5d0*cn(1)*(cn(3)-cn(2))**2/
     $  (cn(3)-cn(2))
      endif
      return
      end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
      subroutine orthrobreg(y,x,yy,res,n,ks,beta,beta0,cn,cpp,sig,ssrho
     $     ,scale,red,cnr)
      integer n,ks
      double precision y(n),x(n,ks),yy(n),res(n),beta(ks),beta0(ks)
      double precision cn,sig,ssrho(3),cpp,cnr(3)
      double precision sr1,sr2,ss,syx,sig0
      double precision rhoh,psih,psih1,rrhoh,rpsih,rpsih1
      integer ic
      logical scale,red
C
C
C
      do 1 i=1,n
         res(i)=y(i)
 1    continue
      sr1=1d10
      ic=0
      sig0=sig
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
         if(red) then
            sr2=sr2+rrhoh(res(i)/sig,cnr)
         else
            sr2=sr2+rhoh(res(i)/sig,cn)
         endif
 50   continue
      if(sr1-sr2.gt.1d-3*sr2) then
         sr1=sr2
         goto 100
      endif
      if(.not.scale) goto 200
      sig=0d0
      do 55 i=1,n
         if(red) then
            sig=sig+rpsih(res(i)/sig0,cnr)**2
         else
            sig=sig+psih(res(i)/sig0,cn)**2
         endif
 55   continue
       sig=dsqrt(sig/(cpp*dble(n-ks)))*sig0
       if(dabs(sig/sig0-1d0).gt.1d-2) then
         sig0=sig
         goto100
      endif
 200  continue

      ssrho(1)=0d0
      ssrho(2)=0d0
      ssrho(3)=0d0
      do 60 i=1,n
         if(red) then
            ssrho(1)=ssrho(1)+rrhoh(res(i)/sig,cnr)
            ssrho(2)=ssrho(2)+rpsih(res(i)/sig,cnr)**2
            ssrho(3)=ssrho(3)+rpsih1(res(i)/sig,cnr)
         else
            ssrho(1)=ssrho(1)+rhoh(res(i)/sig,cn)
            ssrho(2)=ssrho(2)+psih(res(i)/sig,cn)**2
            ssrho(3)=ssrho(3)+psih1(res(i)/sig,cn)
         endif
 60   continue
      return
      end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      subroutine lmmdch(y,x,n,k,xx,xxx,y1,y2,d,r,beta,xinv,ia,intercept
     $     ,ss,nv,ssr,alpha,q)
      integer n,k,q
      double precision y(n),x(n,k),xx(n*k),xxx(n*k),y1(n),y2(n) ,d(k)
     $     ,r(k),beta(k),xinv(k**2),ss(2**k+2),ssr(2**k+2)
      double precision alpha
      integer ia(k),nv(2**k,2)
      logical intercept
C
      double precision ss0,ss1,pval,util1,util2,pv1,pval1
      double precision betai,mn
      integer id,ks,ns,np,ni,qks
      logical inv
C
      inv=.false.
      id=1
      mn=0d0
      ss(1)=0d0
      do 300 i=1,n
         ss(1)=ss(1)+y(i)**2
         mn=mn+y(i)
 300  continue
      mn=mn/dble(n)
      if(intercept) ss(1)=ss(1)-dble(n)*mn**2
C
      if(intercept) then
         kk=2**(k-1)
      else
         kk=2**k
      endif
      ni=0
      if(intercept) ia(k)=1
C
      ni=0
      do 40 iv=1,kk-1
         if(intercept) then
            call decode(iv,k-1,ia)
            ia(k)=1
         else
           call decode(iv,k,ia)
        endif
         ss(iv+1)=0d0       
         ks=0
         do 20 i=1,k
            ks=ks+ia(i)
 20     continue
        call xsubset1(x,xx,n,k,ks,ia,id)
        call lsq(xx,y,xxx,y1,n,ks,d,r,beta,xinv,y2,inv)
        do 31 i=1,n
           ss(iv+1)=ss(iv+1)+y2(i)**2
 31     continue
 40   continue
C

      do 80 iv=1,kk-1
         ns=0
         np=0
         if(intercept) then
            call decode(iv,k-1,ia)
c            ia(k)=1
             do 500 iu=1,k-1
               ns=ns+ia(iu)*2**(iu-1)
               np=np+ia(iu)
 500        continue
          else
            call decode(iv,k,ia)
             do 600 iu=1,k
               ns=ns+ia(iu)*2**(iu-1)
               np=np+ia(iu)
 600        continue
         endif
         ss1=ss(ns+1)
         nst=ns
         ks=0
         pv1=0d0
         do 45 i=1,k
            ks=ks+ia(i)
 45      continue
C
         k1=k
         if(intercept)  k1=k-1
C
        do 70 is=1,k1
             if(ia(is).eq.1) then
               ia(is)=0
               ns=0
               do 50 iu=1,k1
                  ns=ns+ia(iu)*2**(iu-1)
 50            continue
               ss0=ss(ns+1)
               ia(is)=1
               util1=ss1/ss0
               util2=dble(n-ks)/2d0
               if(util1.le.1d-20) then
                  pval1=0d0
               else
                  pval1=betai(util1,util2,0.5d0)
               endif
               qks=q-ks
               pval=betai(pval1,1d0,dble(qks)+1d0)
               if(pval.gt.alpha) goto 80
            endif
 70      continue
         ni=ni+1
         nv(ni,1)=nst
         nv(ni,2)=np
         ssr(ni)=ss(nst+1)
 80   continue
      return
C
      end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
        subroutine roblmmdch(y,x,n,k,alpha,x1,x2,x3,y1,d,r,beta,xinv
     $       ,res,beta0,cn,sig,ssrho,cpp,ia,ib,pp,intercept,nv,ssr,q,
     $       res1,red,cnr)
      integer n,k,q
      double precision y(n),x(n,k),x1(n,k),x2(n*k),x3(n*k),y1(n),d(k)
     $     ,r(k),beta(k),xinv(k**2),res(n),beta0(k),pp(k,2),ssr(2**k),
     $     res1(n),cnr(3)
      double precision alpha,cn,cpp,sig,ssrho(3)
      integer ia(k),ib(k),nv(2**(k+1),2)
      logical intercept,scale,red
C
      integer id,ks,ni,ivv,np
      double precision sig0
C
C
      sig0=sig
      id=1
      ni=0
      if(intercept) then
         kk=2**(k-1)
         q=k-1
      else
         kk=2**k
         q=k
      endif
      do 40 iv=1,kk-1 
         ssr(iv)=0d0
       if(intercept) then
            call decode(iv,k-1,ia)
            ia(k)=1
         else
            call decode(iv,k,ia)
         endif
         ks=0
         do 20 i=1,k
            ks=ks+ia(i)
 20     continue
        if(k.eq.2..and.ks.eq.1) goto 40
        call xsubset1(x,x1,n,k,ks,ia,id)
        sig=sig0
        scale=.true.
        call robregp(y,x1,y1,x2,x3,xinv,n,ks,d,r,beta,res,beta0,cn,sig
     $       ,ssrho,cpp,ib,pp,q+1,intercept,res1,scale,red,cnr)
        kss=ks
        if(intercept) kss=ks-1
       do 30 j=1,kss
           if(pp(j,1).gt.alpha) goto 40
 30     continue
        ni=ni+1
        ivv=0
        np=0
        if(intercept) then
            do 35 ij=1,k-1 
           ivv=ivv+ia(ij)*2**(ij-1)
           np=np+ia(ij)
 35     continue
        else
           do 24 ij=1,k 
           ivv=ivv+ia(ij)*2**(ij-1)
           np=np+ia(ij)
 24     continue
      endif
        nv(ni,1)=ivv
        nv(ni,2)=np
        ssr(ni)=sig
 40   continue 
C
      end
C
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
C
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
c               xx(i+(ic-1)*n)=x(i,j)
 10         continue
         endif
 20   continue
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      subroutine xindsub(x,xx,n,k,ks,ind)
      integer n,k,ks
      double precision x(n,k),xx(n,ks)
      integer ind(ks)
C
      integer ic
C
      ic=0
      do 20 j=1,ks
            do 10 i=1,n
               xx(i,j)=x(i,ind(j))
 10         continue
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      subroutine triggen(n,m,tr)
      integer n,m 
      double precision tr(n,2*m)
C
      double precision pi
C
      pi=4*datan(1d0)
C
      do 10 i=1,n
        tr(i,1)=dsin(pi*dble(i)/dble(n))
        tr(i,2)=dcos(pi*dble(i)/dble(n))
 10   continue
C
      if(m.eq.1) return
C
      jj=3
 20   continue
      do 30 i=1,n
      tr(i,jj)=tr(i,jj-2)*tr(i,2)+tr(i,jj-1)*tr(i,1)
 30   continue
      do 40 i=1,n
      tr(i,jj+1)=tr(i,jj-1)*tr(i,2)-tr(i,jj-2)*tr(i,1)
 40   continue
      if(jj+1.eq.2*m) return
      jj=jj+2
      goto 20
      return
      end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
C
C
      subroutine add2(a,k)
      integer k
      integer a(k)
C
      integer i,is,j,jj,ns
C
C
      is=0
      do 10 i=1,k
         is=is+a(i)
 10   continue
      if(is.eq.k)  then
         do 20 i=1,k
            a(i)=0
 20      continue
         return
      endif
C
      if(a(k).eq.0) then
         j=k-1
         jj=0
 1       continue
         if(a(j).eq.1) then
            jj=j
         elseif(j.ge.2) then
            j=j-1
            goto 1
         endif
         a(jj)=0
         a(jj+1)=1
         return
      endif
C
      j=0
      ns=0
      is=0
      i=k
 2    continue
      if(is.eq.0.and.a(i).eq.1) then
         ns=ns+1
      else
         is=1
      endif
      if(a(i).eq.1.and.is.eq.1) then
         j=i
         goto 3
      elseif(i.ge.2) then
         i=i-1
         goto 2
      endif
C
 3    continue
      if(j.ge.1) then
         do 4 ii=j,k
            a(ii)=0
 4       continue
         do 5 ii=j+1,j+1+ns
            a(ii)=1
 5       continue
      else
         do 6 ii=1,k
            a(ii)=0
 6       continue
         do 7 ii=1,ns+1
            a(ii)=1
 7       continue
         return
      endif
      end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      subroutine allprx(y,x,n,k,xx,xxx,y1,y2,d,r,beta,xinv,ia,intercept,
     $     ss,nv,ssr,alpha,q,kmxx,kmx,ib)
      integer n,k,q,kmx,kmxx
      double precision y(n),x(n,k),xx(n*k),xxx(n*k),y1(n),y2(n) ,d(k)
     $     ,r(k),beta(k),xinv(k**2),ss(kmxx+2),ssr(kmxx+2)
      double precision alpha
      integer ia(k),ib(k),nv(kmxx,2)
      logical intercept
C
      double precision ss0,ss1,pval,util1,util2,pval1
      double precision betai,mn
      integer id,ks,np,ni,nst,kk
      logical inv
C

      inv=.false.
      id=1
      mn=0d0
      ss(1)=0d0
      do 300 i=1,n
         ss(1)=ss(1)+y(i)**2
         mn=mn+y(i)
 300  continue
      mn=mn/dble(n)
      if(intercept) ss(1)=ss(1)-dble(n)*mn**2
C
C
      kk=k
      if(intercept) kk=k-1
C
      do 1 i=1,kk
         ia(i)=0
 1    continue
      if(intercept) ia(k)=1
C
C
      ns=1
 10   continue
      call add2(ia,kk)
      ks=0
      do 20 i=1,kk
         ks=ks+ia(i)
 20   continue
      ks0=ks
      if(ks.eq.0) goto 45
      if(kmx.gt.0.and.ks.gt.kmx) goto 45
      if(intercept) ks=ks+1
      call xsubset1(x,xx,n,k,ks,ia,id)
      call lsq(xx,y,xxx,y1,n,ks,d,r,beta,xinv,y2,inv)
      ns=ns+1
      ss(ns)=0d0
      do 40 i=1,n
         ss(ns)=ss(ns)+y2(i)**2
 40   continue
      goto 10
C
 45   continue
C
      do 2 i=1,kk
         ia(i)=0
 2    continue
      if(intercept) ia(k)=1
C
      ni=0
 50   continue
         call add2(ia,kk)
         ks=0
      do 60 i=1,kk
         ks=ks+ia(i)
 60   continue
      if(ks.eq.0) goto 500
      if(kmx.gt.0.and.ks.gt.kmx) goto 500
      np=ks
      if(intercept) ks=ks+1
      call retn(ia,ib,kk,nst)
      ss1=ss(nst)
C
      do 100 is=1,kk
         if(ia(is).eq.0) goto 100
         ia(is)=0
         call  retn(ia,ib,kk,ns)  
         ss0=ss(ns)
         ia(is)=1
         util1=ss1/ss0
         util1=dmin1(util1,1d0-1d-12)
         util2=dble(n-ks)/2d0
         if(util1.le.1d-20) then
            pval1=0d0
         else
            pval1=betai(util1,util2,0.5d0)
            pval=betai(pval1,1d0,dble(q-ks)+1d0)
            if(pval.gt.alpha)  goto 50
         endif
 100  continue
      ni=ni+1
      nv(ni,1)=nst
      nv(ni,2)=np
      ssr(ni)=ss1
      goto 50
 500  continue
      return
      end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      subroutine retn(ia,ib,k,ns)
      integer k,ns
      integer ia(k),ib(k)
C
      integer id, kj,ns1
C
C
      kj=0
      do 10 i=1,k
         kj=kj+ia(i)
 10   continue
      if(kj.eq.0) then
         ns=1
         return
      endif
C
C
      ns=1
      ns1=1
      do 30 i=1,kj-1
         ns1=(ns1*(k-i+1))/i
         ns=ns+ns1
 30   continue
      ns=ns+1
      do 20 i=1,k
         if(i.le.kj) then
            ib(i)=1
         else
            ib(i)=0
         endif
 20   continue
C
 40   continue
      id=0
      do 50 j=1,k
         if(ia(j).ne.ib(j)) id=id+1
 50   continue
      if(id.gt.0) then
         call add2(ib,k)
         ns=ns+1
         goto 40
      endif
      return
      end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      subroutine retia(ns,ia,k)
      integer ns,k
      integer ia(k)
C
      integer ii,ns1,ns2,ns3
      if(ns.gt.2**k) ns=2**k
      do 10 i=1,k
         ia(i)=0
 10   continue
C
      if(ns.eq.1) return      
      ns2=1
      ns1=1
      ns3=1
      ii=1
 20   continue
      ns2=(ns2*(k-ii+1))/ii
      ns1=ns1+ns2
      if(ns1.le.ns-1) then
         ns3=ns1+1
         ii=ii+1
         goto 20
      endif
C
      do 30 i=1,k
         if(i.le.ii) then
            ia(i)=1
         else
            ia(i)=0
         endif
 30   continue
      if(ns.eq.ns3+1) return
       do 40 i=ns3+2,ns+1
         call add2(ia,k)
 40   continue
      return
      end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      subroutine lagg(x,n,k,lag,xl,y)
      integer n,k,lag
      double precision x(n,k),xl(n-lag,k*lag),y(n-lag)
C
C
      do 10 i=1,n-lag
         y(i)=x(lag+i,1)
 10   continue
C
      do 40 j=1,k
         do 30 lg=1,lag
            do 20 i=1,n-lag
               xl(i,(j-1)*lag+lg)=x(lag-lg+i,j)
 20         continue
 30      continue
 40   continue
      return
      end
