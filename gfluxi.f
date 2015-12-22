      SUBROUTINE GFLUXI(TEMP,TAU,W0,COSBAR,wavenum
     &     ,RSF,FUP,Fdown)
      implicit none

c  mean ubar
      DOUBLE PRECISION UBARI
      INTEGER nlayer, nlevel, ngaussi
      PARAMETER (ngaussi=8, nlayer = 64, nlevel = 65)
      PARAMETER (UBARI = 0.5d0)
c  THIS SUBROUTINE TAKES THE OPTICAL CONSTANTS AND BOUNDARY CONDITONS
C  FOR THE INFRARED FLUX AT ONE WAVELENGTH AND SOLVES FOR THE FLUXES AT
C  THE LEVELS. THIS VERSION IS SET UP TO WORK WITH LAYER OPTICAL DEPTHS
C  MEASURED FROM THE TOP OF EACH LAEYER.     THE TOP OF EACH LAYER HAS  
C  OPTICAL DEPTH ZERO.  IN THIS SUB LEVEL N IS ABOVE LAYER N. THAT IS LAYER N
C  HAS LEVEL N ON TOP AND LEVEL N+1 ON BOTTOM. OPTICAL DEPTH INCREASES
C  FROM TOP TO BOTTOM. SEE C.P. MCKAY, TGM NOTES.
      DOUBLE PRECISION TEMP(nlevel),W0(nlayer),COSBAR(nlayer)
      DOUBLE PRECISION FDOWN(nlevel), FUP(nlevel)
      DOUBLE PRECISION B0(Nlayer),B1(Nlayer),ALPHA(Nlayer),LAMDA(nlayer)
      DOUBLE PRECISION XK1(nlayer),XK2(nlayer)
      DOUBLE PRECISION GAMA(nlayer),CP(nlayer),CM(nlayer),CPM1(nlayer)
      DOUBLE PRECISION CMM1(nlayer),E1(nlayer),E2(nlayer),E3(nlayer)
      DOUBLE PRECISION E4(nlayer)
      DOUBLE PRECISION g(nlayer),xj(nlayer),h(nlayer),xk(nlayer)
      DOUBLE PRECISION alpha1(nlayer)
      DOUBLE PRECISION alpha2(nlayer)
      DOUBLE PRECISION sigma1(nlayer),sigma2(nlayer)
      DOUBLE PRECISION w(ngaussi),x(ngaussi)
      DOUBLE PRECISION fpt(nlevel),fmt(nlevel)
      DOUBLE PRECISION em(nlevel),em2(nlevel),em3(nlevel)
      DOUBLE PRECISION epp(nlevel)
      DOUBLE PRECISION term, rsf,alphax,bottom,bsurf,btop,emm
      DOUBLE PRECISION EP, obj, pi, tautop, ugauss
      DOUBLE PRECISION wavenum
      DOUBLE PRECISION TAU(nlayer)
      INTEGER iflag,j,ng 
      double precision BBPLK
      EXTERNAL BBPLK

c     GAUSS ANGLES AND GAUSS WEIGHTS FOR GAUSSIAN INTEGRATION
c     MOMENTS (USE FIRST MOMENT VALUES) N=3
C
c      x = [0.2123405382d0, 0.5905331356d0,0.9114120405d0]
c      w = [0.0698269799d0, 0.2292411064d0, 0.2009319137d0]
C
C     GAUSS ANGLES AND WEIGHTS FOR GAUSSIAN INTEGRATION MOMENTS
C     (USE FIRST MOMENT ONLY)  N=8
C     X = gauss angle
        x =  [0.0446339553d0, 0.1443662570d0,
     &     0.2868247571d0, 0.4548133152d0, 0.6280678354d0,
     &     0.7856915206d0, 0.9086763921d0, 0.9822200849d0]
C     w = gauss weight
      w =  [0.0032951914d0, 0.0178429027d0,
     &     0.0454393195d0, 0.0791995995d0, 0.1060473594d0,
     &     0.1125057995d0, 0.0911190236d0, 0.0445508044d0]
c       5 Gauss points for mu integration:
c     DATA GANGLE /0.0985350858, 0.3045357266, 0.5620251898, 0.8019865821,
c    &           0.9601901429/
c     DATA GWEIGHT/0.0157479145, 0.0739088701, 0.1463869871,  0.1671746381,
c    &           0.0967815902/
      iflag = 0
      PI = 3.14159265358979323846d0
      do  j =1, nlevel
         fdown(j) = 0.0d0
         fup(j) = 0.0d0
      end do
      
      DO 20 J=1,NLAYER
         ALPHA(J)=DSQRT( (1.d0-W0(J))/(1.d0-W0(J)*COSBAR(J)))
         LAMDA(J)=ALPHA(J)*(1.d0-W0(J)*COSBAR(J))/UBARI
         GAMA(J)=(1.d0-ALPHA(J))/(1.d0+ALPHA(J))
         term=0.5d0/(1.d0-w0(j)*cosbar(j))
         B1(J)=(BBPLK(wavenum,TEMP(J+1)) - 
     &        BBPLK(wavenum,TEMP(J))) /(tau(j))
         B0(J)=BBPLK(wavenum,TEMP(J))
         if (tau(j).lt.1e-6) then
            b1(j) = 0.0
            b0(j) = 0.5d0*(BBPLK(wavenum,temp(j))
     &           +BBPLK(wavenum,temp(j+1)))
         endif
         cp(j)=b0(j)+b1(j)*tau(j)+b1(j)*term
         cm(j)=b0(j)+b1(j)*tau(j)-b1(j)*term
         cpm1(j)=b0(j)+b1(j)*term
         cmm1(j)=b0(j)-b1(j)*term

 20   CONTINUE
      
C     * NOW CALCULATE THE EXPONENTIAL TERMS NEEDED
C     * FOR THE TRIDIAGONAL ROTATED LAYERED METHOD
C     * WARNING IF tau(j) IS GREATER THAN ABOUT 35
C     * WE CLIP IT TO AVOID OVERFLOW.
C     * EXP (TAU) - EXP(-TAU) WILL BE NONSENSE THIS IS
C     * CORRECTED IN THE DSOLVER ROUTINE. ?FLAG? 
      DO 8 J=1,NLAYER
         EP=EXP(35.d0)
         IF (LAMDA(J)*tau(j) .LT. 35.)  EP=EXP(LAMDA(J)*tau(j))
         EMM=1.d0/EP
c     This is the method used by McKay, it does not rotate the e's
         E1(J)=EP+GAMA(J)*EMM
         E2(J)=EP-GAMA(J)*EMM
         E3(J)=GAMA(J)*EP+EMM
         E4(J)=GAMA(J)*EP-EMM
c     write (*,9999) j,e1(j),e2(j),e3(j),e4(j)
c     This is the method actually printed in Toon et al.  It does rotate the e's
c     E1(J)=1.d0+GAMA(J)*EMM
c     E2(J)=1.d0-GAMA(J)*EMM
c     E3(J)=GAMA(J)+EMM
c     E4(J)=GAMA(J)-EMM
c     if (j.eq.25) print *,'es',e1(j),e2(j),e3(j),e4(j)
 8    CONTINUE
      
      TAUTOP=TAU(1)
      BTOP=(1.d0 - EXP(-TAUTOP/ubari))*BBPLK(wavenum,TEMP(1))
      BSURF = BBPLK( wavenum, TEMP(nlevel)) 
      bottom = (bsurf+b1(nlayer)*ubari)
c     print *,'bsurf,bottom',bsurf,bottom,pi,b1(nlayer),ubari
      
      CALL DSOLVER(NLAYER,GAMA,CP,CM,CPM1,CMM1
     &     ,E1,E2,E3,E4,BTOP,bottom,RSF,XK1,XK2)
      

     
C     Loop over the ugauss beginning here
      DO 315 NG = 1,ngaussi
        ugauss = x(ng)
c     print *,'mu info',ng,ugauss,w(ng)
        do j=1,nlayer
           if (w0(j).ge.0.01d0) then
c     McKay
              alphax = ((1.d0-w0(j))/(1.d0-w0(j)*cosbar(j)))**0.5d0
              g(j)=2.d0*pi*w0(j)*xk1(j)*(1.d0+cosbar(j)*alphax)/
     &           (1.d0+alphax)
              h(j)=2.d0*pi*w0(j)*xk2(j)*(1.d0-cosbar(j)*alphax)/
     &             (1.d0+alphax)
              xj(j)=2.d0*pi*w0(j)*xk1(j)*(1.d0-cosbar(j)*alphax)/
     &             (1.d0+alphax)
              xk(j)=2.d0*pi*w0(j)*xk2(j)*(1.d0+cosbar(j)*alphax)/
     &           (1.d0+alphax)
              alpha1(j)=2.d0*pi*(b0(j)
     &           + b1(j)*(ubari*w0(j)*cosbar(j)/(1.d0-w0(j)*cosbar(j))))
              alpha2(j)=2.d0*pi*b1(j)
              sigma1(j)=2.d0*pi*(b0(j)
     &           - b1(j)*(ubari*w0(j)*cosbar(j)/(1.d0-w0(j)*cosbar(j))))
              sigma2(j)=alpha2(j)
           else 
              g(j)=0.0
              h(j)=0.0
              xj(j)=0.0
              xk(j)=0.0
              alpha1(j)=2.d0*pi*b0(j)
              alpha2(j)=2.d0*pi*b1(j)
              sigma1(j)=alpha1(j)
              sigma2(j)=alpha2(j)
           endif
        enddo
        
        fpt(nlevel)=2.d0*pi*(bsurf+b1(nlayer)*ugauss)
        fmt(1)=2.d0*pi*(1.d0 - EXP(-TAUTOP/ugauss))*
     &       BBPLK(wavenum,TEMP(1))
        do j=1,nlayer
           em(j)=EXP(-LAMDA(J)*tau(j))
           em2(j)=EXP(-tau(j)/ugauss)
           em3(j)=em(j)*em2(j)
           epp(j)=exp(35.d0)
           obj=lamda(j)*tau(j)
           if (obj .lt. 35.) epp(j)=exp(obj)
           
           fmt(j+1)=fmt(j)*em2(j) + 
     %          (xj(j)/(lamda(j)*ugauss+1.d0))*(epp(j)-em2(j)) +
     %          (xk(j)/(lamda(j)*ugauss-1.d0))*(em2(j)-em(j)) +
     %          sigma1(j)*(1.d0-em2(j)) +
     %          sigma2(j)*(ugauss*em2(j)+tau(j)-ugauss)
        enddo
        do j=nlayer,1,-1
           fpt(j)=fpt(j+1)*em2(j) + 
     %          (g(j)/(lamda(j)*ugauss-1.d0))*(epp(j)*em2(j)-1.d0) +
     %          (h(j)/(lamda(j)*ugauss+1.d0))*(1.d0-em3(j)) +
     %          alpha1(j)*(1.d0-em2(j)) +
     %          alpha2(j)*(ugauss-(tau(j)+ugauss)*em2(j))
        enddo

        
c     print *,' angle ratios',ng,x(ng),w(ng),gweight 
c        DO 60 J=1,nlevel
c     gauss angle integration goes from 0 to 1
c     azimuthal intgration (in fpt already) accounts
c     for the other side
           FUP(1)=FUP(1)+w(ng)*fpt(1)
c     FDOWN(J)=FDOWN(J)+w(ng)*fmt(j)
c 60     CONTINUE
        
        
 315  continue

      RETURN
      END

      
      double precision function bbplk( waven, T )
      implicit none
      
      double precision T, waven, wavelen      
      double precision C, h, pi, kb

      parameter (c = 299792458.d0 , kb = 1.38064852d-23,
     &     h = 6.62607004d-34, pi =  3.14159274d0)
      
      wavelen = 1.d-6 * (1.0d4 / waven)
      
      
      bbplk = 1.d-6 * ((2.d0*h * c**2) / wavelen**5) /
     &     (exp(h*c/(wavelen*kb*T)) - 1.d0)
      
      return
      end     

