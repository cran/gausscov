C
C
C
      subroutine fstepwise(y,x,n,k,x2,res,ia,alpha,kmx,pp,kex
     $     ,intercept,minss,ss01,qq,kmn)
      integer n,k,kmn,qq,kmx
      double precision y(n),x(n,k),x2(n),res(n) ,pp(k+1,2),minss(k),
     $     ss01(k)
      integer ia(k),kex(k)
      logical intercept
      double precision alpha
C
      integer icount,ks,ic,ik,kr,nex
      double precision ss0,ss1,amss1,pval,util1,util2,nx1,b,cf,pval1,mi
     $     ,ssy,ssx
      double precision betai
C
      if(intercept)  kmx=kmx+1
      ic=0 
      do 1 j=1,k
         ia(j)=0
 1    continue
C
      nex=0
      do 11 ik=1,k
         if(kex(ik).gt.0) then
            ia(kex(ik))=1
            nex=nex+1
         endif
 11   continue
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
      kr=0
      do 19 ik=1,k
         if(ia(ik).eq.1) kr=kr+1
 19    continue
      if(qq.eq.0) then
         kr=k-kr
      else
         kr=qq-kr
      endif
C
 2    continue

      if(ks.eq.k) goto 600
C
C
C
      ks=ks+1
      amss1=ss0
      do 40 kk=1,k
         if(ia(kk).eq.1)  goto 40
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
      if(amss1.lt.1d-10) then
         pval=0d0
         icount=icount+1
         pp(icount,1)=dble(ic)
         pp(icount,2)=0d0
         minss(icount)=0d0
         ss01(icount)=0d0
         kmx=icount
         return
      else
         util1=amss1/ss0
         util2=dble(n-icount-1)/2d0   
         pval1=betai(util1,util2,0.5d0)
         pval=betai(pval1,1d0,dble(kr+2-icount)-1d0)
      endif
      if(kmx.gt.0.and.icount.ge.kmx) goto 600
      if(icount.lt.kmn) goto 49
      if(pval.gt.alpha.and.icount.ge.kmn) goto 600
 49   continue
      icount=icount+1
      pp(icount,1)=dble(ic)
      pp(icount,2)=pval
      minss(icount)=amss1
      ss01(icount)=amss1/ss0
      if(kmx.gt.0.and.icount.ge.kmx)  then
         kmx=icount
         goto 600
      endif
      ia(ic)=1
      b=0d0
      nx1=0d0
      do 45 i=1,n
         b=b+x2(i)*res(i)
         nx1=nx1+x2(i)**2
 45   continue
      b=b/nx1
      ss0=0d0
      cf=dsqrt(dble(n)/nx1)
      nx1=0d0
      do 50 i=1,n
         res(i)=res(i)-b*x2(i)
         ss0=ss0+res(i)**2
         x2(i)=cf*x2(i)
         nx1=nx1+x2(i)**2
 50   continue
      if(icount+nex.eq.k) then
         kmx=icount
         return
      endif
      do 60 kk=1,k
         if(ia(kk).eq.1) goto 60
         b=0d0
         do 61 i=1,n
            b=b+x(i,kk)*x2(i)
 61      continue
         b=b/dble(n)
         ssx=0d0
         do 62 i=1,n
            x(i,kk)=x(i,kk)-b*x2(i)
            ssx=ssx+x(i,kk)**2
 62      continue
         if(ssx.lt.1d-10) ia(kk)=1
 60   continue
      goto 2
 600  continue
      kmx=icount
      return
      end
C
C
C    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
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
c      subroutine edge(edg,ne,kmax,ke,ned)
c      integer ne,kmax,ned
c      integer edg(ne,2),ke(kmax)
C
c      integer i0,i1,i,j,ll
C
C
c      do 1 i=1,ne
c         if(edg(i,1).gt.edg(i,2)) then
c            ll=edg(i,1)
c            edg(i,1)=edg(i,2)
c            edg(i,2)=ll
c         endif
c 1      continue
         
c      call iquicksort(edg,ne,2,1)
c      return
C
c      i0=1
c      i1=1
c      i2=1
c      i=1
c      j=edg(1,1)
c      ic=0
c 10   continue
c      if(edg(i,1).eq.j) then
c         ic=ic+1
c         ke(ic)=edg(i,2)
c         if(i.eq.ne) goto 15
c         if(i.lt.ne) then
c            i=i+1
c            goto 10
c         endif
c      endif
c 15   continue
c      if(ic.gt.kmax) return
c      call iquicksort(ke,ic,1,1)
c      edg(i2,2)=ke(1)
c      edg(i2,1)=j
c      do 20 ii=2,ic
c         if(ke(ii).eq.ke(ii-1)) goto 20
c         i2=i2+1
c         edg(i2,2)=ke(ii)
c         edg(i2,1)=edg(i0,1)
c 20   continue
c      if(i.eq.ne) then
c         ned=i2
c         return
c      endif
c      i2=i2+1
c      if(i.lt.ne) then
c         i0=i
c         j=edg(i,1)
c         ic=0
c         goto 10
c      endif
c      return
c      end
C
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C
C
C
      subroutine graphst(xxx,x,n,k,y,x2,res,ia,alpha,kmx,pp,grph
     $     ,ne,kexc,xinr,minss,nedge,ss01,kmn,lin,iind,grphp)
      integer n,k,kmn,ne,nedge,kmx,lin
      double precision xxx(n,k),x(n,k),y(n),x2(n),res(n),pp(k+1,2)
     $     ,minss(k),ss01(k),grphp(nedge)
      integer ia(k),grph(nedge,3),kexc(k),iind(lin)
      logical xinr
      double precision alpha
C
      integer qq,kmx1,jj

C
C
      qq=k
      ne=0
      kk=k
      if(xinr) kk=k-1
      do 6 il=1,lin
         jj=iind(il)
         do 4 j1=1,k
            do 3 ii=1,n
                  x(ii,j1)=xxx(ii,j1)
 3          continue
 4       continue
         do 10 i=1,n
            y(i)=xxx(i,jj)
 10      continue
         kmx1=kmx
         kexc(1)=jj
         call fstepwise(y,x,n,k,x2,res,ia,alph a,kmx1,pp,kexc
     $        ,xinr,minss,ss01,qq,kmn)
         if(kmx1.eq.0) goto 6
         if(kmx1.eq.1.and.idnint(pp(1,1)).eq.0) goto 6
         do 15, ij=2,kmx1
           if(idnint(pp(ij,1)).ge.1) then
              ne=ne+1
              grph(ne,1)=jj
             grph(ne,2)=idnint(pp(ij,1))
             grphp(ne)=pp(ij,2)
               if(ne.ge.nedge) return
            endif
 15      continue
 6    continue
C
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
      subroutine graphstst(xxx,x,n,k,y,x2,res,ia,alpha,kmx,pp,grph,
     $     ne,kexc,xinr,minss,nedge,ss01,rgrph,kmn,lin,iind)
      integer n,k,kmn,ne,nedge,kmx,lin
      double precision xxx(n,k),x(n,k),y(n),x2(n),res(n),pp(k+1,2)
     $     ,minss(k),ss01(k),rgrph(nedge)
      integer ia(k),grph(nedge,3),kexc(k),iind(lin)
      double precision alpha
C
      logical xinr
      integer ij,kmx1,qq,ek,iq,la,ij0

C
      qq=0
      ne=0
      ek=0
      kk=k
      if(xinr) then
         kk=k-1
      endif
      do 100 il=1,lin
         jj=iind(il)
         ek=0
         do 2 i=1,n
            y(i)=xxx(i,jj)
 2       continue
         do 6 i=1,k
            ia(i)=0
            kexc(i)=0
 6       continue
         kexc(1)=jj
         ek=1
         la=0
 7       continue
         iq=0
         do 8, i=1,k
            if(kexc(i).gt.0) iq=iq+1
 8       continue
         if(iq.eq.kk) goto 100
         do 5 j1=1,k
            do 4 ii=1,n
                  x(ii,j1)=xxx(ii,j1)
 4          continue
 5       continue
         kmx1=kmx
          call fstepwise(y,x,n,k,x2,res,ia,alpha,kmx1,pp,kexc
     $        ,xinr,minss,ss01,qq,kmn)
         if(kmx1.le.0) goto 100
         if(kmx1.eq.1.and.xinr) goto 90
         la=la+1
         ij0=1
         if(xinr) ij0=2
         do 15, ij=ij0,kmx1
            if(xinr.and.idnint(pp(ij,1)).eq.k) goto 15
            if(idnint(pp(ij,1)).ge.1) then
               ne=ne+1
               ek=ek+1
               kexc(ek)=idnint(pp(ij,1))
               grph(ne,1)=jj
               grph(ne,2)=la
               grph(ne,3)=idnint(pp(ij,1))
               rgrph(ne)=pp(ij,2)
               endif
               if(ne.ge.nedge) return
               if(ne.gt.k*n) return
 15      continue
         goto 7
 90      continue
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
               qks=q-np+1
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
c      subroutine xindsub(x,xx,n,k,ks,ind)
c      integer n,k,ks
c      double precision x(n,k),xx(n,ks)
c      integer ind(ks)
C
c      integer ic
C
c      ic=0
c      do 20 j=1,ks
c            do 10 i=1,n
c               xx(i,j)=x(i,ind(j))
c 10         continue
c 20   continue
c      end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      subroutine genint(x,xx,n,k,kk,kex,ord,ind,ji)
      integer n,k,kk,ord,ji
      integer kex(kk,ord),ind(ord)
      double precision x(n,k),xx(n,kk)
C
C
      do 10 i=1,ord
         ind(i)=1
 10   continue
C        
      ji=0
      do 40 j=1,kk
C
         do 30 i=1,n
            xx(i,j)=1d0
            do 20 jj=1,ord
               xx(i,j)=xx(i,j)*x(i,ind(jj))
 20         continue
 30      continue
         ji=ji+1
         do 37 io=1,ord
            if(ind(io).eq.k) then
               kex(ji,io)=0
            else
               kex(ji,io)=ind(io)
            endif
 37      continue
         call inact(ind,k,ord)
 40   continue
      return
      end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
c      subroutine degenint(ic,m,k,ord,ind,iut)
c      integer m,k,ord
c      integer ind(m,ord),ic(m),iut(ord)
C
C
c      do 20 j=1,m
c         call degenint1(ic(j),k,ord,iut)
c         do 10 i=1,ord
c            ind(j,i)=iut(i)-1
c 10      continue
c 20   continue
c      return
c      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
c      subroutine degenint1(ic0,k,ord,ind)
c      integer ic0,k,ord
c      integer ind(ord)
C
C
c      do 10 i=1,ord
c         ind(i)=1
c 10   continue
C
c      ic=0
c 40   continue
c      ic=ic+1
c      if(ic.eq.ic0) return
c      call inact(ind,k,ord)
c      goto 40
c      return
c      end
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
      subroutine lagg(x,n,k,lag,xl,y,jj)
      integer n,k,lag,jj
      double precision x(n,k),xl(n-lag,k*lag),y(n-lag)
C
C
      do 10 i=1,n-lag
         y(i)=x(lag+i,jj)
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
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C     THIS SUBROUTINE PRODUCES GAUSSIAN N(0,1) RANDOM VARIABLES
C
      subroutine gaussrnd(x,n,idum)
      integer n,idum
      double precision x(n)
C
      DOUBLE PRECISION PI,RR,THETA,z(1)
C
      parameter (PI=3.141592654D0)
c      idum=1
C
      do 10 i=1,n,2
         call runif2(1,z,idum)
         rr=dsqrt(-2.0d0*dlog(z(1)))
         call runif2(1,z,idum)
         theta=2.0d0*pi*z(1)
         x(i)=rr*dsin(theta)
         if (i+1.le.n) x(i+1)=rr*dcos(theta)
 10   continue
      return
      end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      subroutine runif2(n,x,idum)
      integer n
      double precision x(n)
      INTEGER IDUM
C
      double precision ran2

C
c      idum=1
      do 10 i=1,n
         x(i)=ran2(idum)
 10   continue
      return
      end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      double precision function ran2(idum)
      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      double precision am,eps,rnmx
      parameter (im1=2147483563,im2=2147483399,am=1d0/dble(im1),
     $     imm1=im1-1,ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,
     $     ir2=3791,ntab=32,ndiv=67108862,eps=1.2d-7,rnmx=1d0-eps)
      integer idum2,j,k,iv(ntab),iy
      save iv,iy,idum2
      data idum2/123456789/,iv/ntab*0/,iy/0/
C
      if(idum.le.0) then
         idum=max0(-idum,1)
         idum2=idum
         do 11 j=ntab+8,1,-1
            k=idum/iq1
            idum=ia1*(idum-k*iq1)-k*ir1
            if(idum.lt.0) idum=idum+im1
            if(j.le.ntab) iv(j)=idum
 11      continue
         iy=iv(1)
      endif
      k=idum/iq1
      idum=ia1*(idum-k*iq1)-k*ir1
      if(idum.lt.0) idum=idum+im1
      k=idum2/iq2
      idum2=ia2*(idum2-k*iq2)-k*ir2
      if(idum2.lt.0) idum2=idum2+im2
      j=1+iy/ndiv
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1) iy=iy+imm1
      ran2=dmin1(am*dble(iy),rnmx)
      return
      end
C
C     ran2 from Numerical Recipes in Fortran 77
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      double precision function rgaus2(idum)
      integer idum
C
      integer iset
      double precision fac,gset,rsq,v1,v2,ran2
      save iset,gset
      data iset/0/
C
      if(idum.lt.0) iset=0
      if(iset.eq.0) then
 1       v1=2d0*ran2(idum)
         v2=2d0*ran2(idum)
         rsq=v1**2+v2**2
         if(rsq.ge.1d0.or.rsq.eq.0d0) goto 1
         fac=dsqrt(-2*dlog(rsq)/rsq)
         gset=v1*fac
         rgaus2=v2*fac
         iset=1
      else
         rgaus2=gset
         iset=0
      endif
      return
      end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C     THIS SUBROUTINE PRODUCES GAUSSIAN N(0,1) RANDOM VARIABLES
C
      subroutine gaussrnd2(x,n,idum)
      integer n,idum
      double precision x(n)
C
      double precision rgaus2
C
      do 10 i=1,n
         x(i)=rgaus2(idum)
 10   continue
      return
      end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
