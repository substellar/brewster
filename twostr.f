      REAL FUNCTION chapman( lc, nfac, fac, dtauc )
*
*     Calculates the Chapman-factor, Eq. B2 (DS).
*
* I n p u t       v a r i a b l e s:
*
*      lc        : Computational layer
*      dtauc     : Optical thickness of layer lc (un-delta-m-scaled)
*      fac       : Geomtrical correction factor, calculated by routine
*                  geofac.
*      nfac      : Number of terms in Eq. B2 (DS)
*
* O u t p u t      v a r i a b l e s:
*
*      chapman   : Chapman-factor. In a pseudo-spherical atmosphere 
*                  replace EXP( -tau/umu0 ) by EXP( -chapman ) in the
*                  beam source. 
*
* I n t e r n a l     v a r i a b l e s:
*
*      fact      : =1 for first sum in Eq. B2 (DS)
*                  =2 for second sum in Eq. B2 (DS)
*
c +-------------------------------------------------------------------+
c
c  LOCAL SYMBOLIC DIMENSIONS 
c
c  MXCLY  = Max no. of computational layers
c  MXULV  = Max no. of output levels
c  MXCMU  = Max no. of computation polar angles
c  MXUMU  = Max no. of output polar angles
c  MXPHI  = Max no. of output azimuthal angles
c +-------------------------------------------------------------------+

      INTEGER MXCLY, MXULV, MXTCMU, MI9M2, NNLYRI
      INCLUDE "DISORT.MXD"
      PARAMETER ( MXTCMU = 2, MI9M2 = 4,
     &          NNLYRI = MXTCMU*MXCLY)
*     
      INTEGER nfac(*)
      REAL fac(mxcly,*), dtauc(*)
*     
      sum = 0.0
      IF ( nfac(lc) .LT. 0 )  THEN
         sum = 1.E+20
      ELSE
         DO j = 1, nfac(lc)
*     Include factor of 2 for zenang .gt. 90, second sum in Eq. B2 (DS)
            fact = 1.0
            IF( j.GT.lc ) fact = 2.0
            sum = sum + dtauc(j)*fact*fac(lc,j)
         ENDDO
         IF( nfac(lc).GT.lc ) THEN
            sum = sum + dtauc(lc) * fac(lc,j+1)
         END IF
      ENDIF
      chapman = sum
      RETURN
      END
*
      SUBROUTINE geofac( fac, nfac, z_lay, nlyr, zd, zenang, r )
*
* Calculates the geometric correction factor needed for spherical
* geometry. See Appendix B, (DS).
*
* I n p u t       v a r i a b l e s:
*
*      lc        : Computational layer
*      nlyr      : Number of layers in atmospheric model
*      zd(lc)    : lc = 0, nlyr. zd(lc) is distance from bottom 
*                  surface to top of layer lc. zd(nlyr) = 0.0 km
*      zenang    : Solar zenith angle as seen from bottom surface
*      r         : Radius of earth. NOTE: Use the same dimension as zd,
*                  for instance both in km.
*      z_lay     : Where in the layer the Chapman function is to be
*                  computed. E.g. 0.0, bottom of layer, 0.5, middle
*                  of layer and 1.0 top of layer.
*
* O u t p u t      v a r i a b l e s:
*
*      fac        : geometrical correction factor
*      nfac       : Number of terms in Eq. B2 (DS)
*
* I n t e r n a l     v a r i a b l e s:
*
*      dhj       : delta-h-sub-j in Eq. B2 (DS)
*      dsj       : delta-s-sub-j in Eq. B2 (DS)
*      fact      : =1 for first sum in Eq. B2 (DS)
*                  =2 for second sum in Eq. B2 (DS)
*      rj        : r-sub-j in Eq. B1 (DS)
*      rjp1      : r-sub-j+1 in Eq. B1 (DS)
*      xpsinz    : The length of the line OG in Fig. 1, (DS)
*
*
c +-------------------------------------------------------------------+
c
c  LOCAL SYMBOLIC DIMENSIONS 
c
c  MXCLY  = Max no. of computational layers
c  MXULV  = Max no. of output levels
c  MXCMU  = Max no. of computation polar angles
c  MXUMU  = Max no. of output polar angles
c  MXPHI  = Max no. of output azimuthal angles
c +-------------------------------------------------------------------+

      INTEGER MXCLY, MXULV, MXTCMU, MI9M2, NNLYRI
      INCLUDE "DISORT.MXD"
      PARAMETER ( MXTCMU = 2, MI9M2 = 4,
     &          NNLYRI = MXTCMU*MXCLY)

      INTEGER nfac(*)
      REAL fac(mxcly,*), zd(0:*)
*
      REAL dsj, dhj
*
      pi     = 2.0 * ASIN( 1.0 )
      zenrad = zenang * pi / 180.0
*
      DO lc = 1, nlyr
         nfac(lc) = 0
      ENDDO
*
      DO lc = 1, nlyr
         xp     = r +  zd(lc) + (zd(lc-1) - zd(lc) ) *
     $        z_lay
         xpsinz = xp * SIN( zenrad )
         IF( (zenang.GT.90.0) .AND. (xpsinz.LT.r) ) THEN
            nfac(lc) = -1            
         ELSE
*
*     Find index of layer in which the screening height lies
*
            id = lc
            IF( zenang.GT.90.0 ) THEN
               DO j = lc, nlyr
                  IF( (xpsinz.LT.( zd(j-1) + r ) ) .AND.
     $                 (xpsinz.GE.( zd(j) + r )) ) id = j
               ENDDO
            END IF
*
            DO j = 1, id
               fact2= 1.0
*     
* Include factor of 2 for zenang .gt. 90, second sum in Eq. B2 (DS)
*
               IF( j.GT.lc ) fact2 = 2.0
               IF(j.EQ.id .AND. id.EQ.lc .AND. zenang.GT.90.0) 
     $              fact2 = -1.0
*     
               rj = r + zd(j-1)
               rjp1 = r + zd(j)
               IF(j.EQ.lc .AND. id.EQ.lc) rjp1 = xp
*     
               dhj = zd(j-1) -zd(j)
               IF(id.GT.lc .AND. j.EQ.id) THEN
                  dsj = SQRT(rj*rj - xpsinz*xpsinz )
               ELSE
                  dsj = SQRT( rj*rj - xpsinz*xpsinz ) -
     $                 fact2 * SQRT( rjp1*rjp1 - xpsinz*xpsinz )
               END IF
*     
               fac(lc,j) = dsj / dhj
*     
            ENDDO
            nfac(lc) = id
*     
* Third term in Eq. B2 (DS)
*
            IF( id.GT.lc ) THEN
               dhj = zd(lc-1) -zd(lc)
               dsj = SQRT( xp*xp - xpsinz*xpsinz ) -
     $              SQRT( (zd(lc)+r)*(zd(lc)+r) - xpsinz*xpsinz )
               fac(lc,j+1) = dsj / dhj
            END IF
*
            ENDIF
*
      ENDDO
*
      RETURN
      END
*
      SUBROUTINE hopsol( biggest, ch, chtau, cmu, fbeam, ggprim, kk,
     $     ncut, oprim, pi, pkag, pkagc, planck, rr,
     $     smallest, taucpr, umu0, xb_0d, xb_0u, xb_1d, xb_1u,
     $     xp_0,xp_1, yb_0d, yb_0u, yb_1d, yb_1u, yp_0d, yp_0u,
     $     yp_1d, yp_1u, zb_a, zp_a )
*
* Calculates the homogenous and particular solutions to the
* radiative transfer equation in the two-stream approximation,
* for each layer in the medium.
*
*    I n p u t     v a r i a b l e s:
*
*       ch       :  Chapman correction factor
*       cmu      :  Abscissae for gauss quadrature over angle cosine
*       ncut     :  Number of computational layer where absorption
*                     optical depth exceeds -abscut-
*       oprim    :  Delta-m scaled single scattering albedo
*       pkag,c   :  Planck function in each layer
*       taucpr   :  Cumulative optical depth (delta-m-scaled)
*       (remainder are 'twostr' input variables)
*
*   O u t p u t     v a r i a b l e s:
*
*       kk       :  Eigenvalues 
*       rr       :  Eigenvectors at polar quadrature angles 
*       xp_0,    :  Thermal source coefficients x-sub-zero and
*        xp_1         x-sub-one in Eq.  KST(22)
*       yb_0d,u, :  Thermal source vectors, Eqs. KST(24-25)
*        yb_1d,u     
*       yp_0d,u, :  Beam source vectors, Eq. KST(24-25)
*        yp_1d,u     
*       zb_a     :  Beam source coefficient alfa in Eq. KST(22)
*       zp_a     :  Thermal source coefficient alfa in Eq. KST(22)
*
*
      LOGICAL planck
      REAL ch(*), chtau(0:*), ggprim(*), kk(*), large, oprim(*),
     $      pkag(0:*), pkagc(*), rr(*), taucpr(0:*), xb_0d(*),
     $     xb_0u(*), xb_1d(*), xb_1u(*), xp_0(*), xp_1(*), yb_0d(*),
     $     yb_0u(*), yb_1d(*), yb_1u(*), yp_0d(*), yp_0u(*), yp_1d(*),
     $     yp_1u(*), zb_a(*), zp_a(*)
*
* The calculation of the particular solutions require some care, small
* and big have been set so that no problems should occurr on 32-bits
* machine running single precision
*
      big   = sqrt(biggest) / 1.E+10
      small = 1.E+30*smallest
CBM  moved the definition of large here - was in the 
CBM  IF ( fbeam.GT.0.0 ) loop before; large was used un-initialized 
CBM  therefore in thermal calculations!
      large = LOG(biggest)-20. 
*
* ===================  Begin loop on computational layers  =============
*
      DO 100  LC = 1, NCUT
         ls = 2*(lc-1) + 1
*
* Calculate eigenvalues -kk- and eigenvector -rr-, Eqs. KST(20-21) 
*
         beta   = 0.5 * ( 1.-3.*ggprim(lc)*cmu*cmu )
         fact1  = 1. - oprim(lc)
         fact2  = 1. - oprim(lc) + 2.*oprim(lc)*beta 
         kk(lc) = (1./cmu) * sqrt( fact1 * fact2 )
         rr(lc) = ( sqrt(fact2)-sqrt(fact1) ) /
     $                                    ( sqrt(fact2)+sqrt(fact1) )
*aky         write(*,*)"eigen", lc, fact2, fact1, beta, ggprim(lc), rr(lc)
*
         IF ( fbeam.GT.0.0 ) THEN
*
* Set coefficients in KST(22) for beam source
*
           q_1 = (fbeam/(4.*pi))*oprim(lc)*(1.-3.*ggprim(lc)*cmu*umu0)
           q_2 = (fbeam/(4.*pi))*oprim(lc)*(1.+3.*ggprim(lc)*cmu*umu0) 
*
           IF ( umu0 .GE. 0.0 ) THEN
              qq = q_2
           ELSE
              qq = q_1
           ENDIF
*
           q0a = EXP(-chtau(lc-1) )
           q0 = q0a*qq
           IF ( q0 .LE. small) THEN
              q2a = 0.0
           ELSE
              q2a = EXP(-chtau(lc) )
           ENDIF
           q2 = q2a*qq
*
* Calculate alpha coefficient 
*
           deltat = taucpr(lc) - taucpr(lc-1)
           zb_a(lc) = 1./ch(lc)
CBM           large = LOG(biggest)-20. 
*
* Dither alpha if it is close to an eigenvalue
*
           denomb =  fact1 * fact2 - (zb_a(lc)*cmu)*(zb_a(lc)*cmu)
           IF ( denomb .LT. 1.E-03 ) THEN
              zb_a(lc) = 1.02*zb_a(lc)
           ENDIF
           q0 = q0a * q_1
           q2 = q2a * q_1
*
* Set constants in Eq. KST(22)
*
           IF ( deltat .LT. 1.E-07 .OR.
     $          ABS(zb_a(lc)*taucpr(lc-1)) .GT. large .OR.
     $          ABS(zb_a(lc)*taucpr(lc)) .GT. large ) THEN
              xb_1d(lc) = 0.0
           ELSE
              xb_1d(lc) = (1./deltat)*(q2*EXP(zb_a(lc)*taucpr(lc)) 
     $             -q0*EXP(zb_a(lc)*taucpr(lc-1)))
           ENDIF
           IF( ABS(zb_a(lc)*taucpr(lc-1)) .GT. large ) THEN
              xb_0d(lc) = q0 
           ELSE
              xb_0d(lc) = q0 * EXP(zb_a(lc)*taucpr(lc-1)) -
     $             xb_1d(lc)*taucpr(lc-1)
           ENDIF
           q0 = q0a * q_2
           q2 = q2a * q_2
           IF ( deltat .LT. 1.E-07.OR.
     $          ABS(zb_a(lc)*taucpr(lc-1)) .GT. large .OR.
     $          ABS(zb_a(lc)*taucpr(lc)) .GT. large ) THEN
              xb_1u(lc) = 0.0
           ELSE
              xb_1u(lc) = (1./deltat)*(q2*EXP(zb_a(lc)*taucpr(lc)) 
     $             -q0*EXP(zb_a(lc)*taucpr(lc-1)))
           ENDIF
           IF( ABS(zb_a(lc)*taucpr(lc-1)) .GT. large ) THEN
              xb_0u(lc) = q0 
           ELSE
              xb_0u(lc) = q0 * EXP(zb_a(lc)*taucpr(lc-1)) -
     $             xb_1u(lc)*taucpr(lc-1)
           ENDIF
*     
* Calculate particular solutions for incident beam source in 
* pseudo-spherical geometry, Eqs. KST(24-25)
*
           denomb =  fact1 * fact2 - (zb_a(lc)*cmu)*(zb_a(lc)*cmu)
           yb_1d(lc) = (  oprim(lc)*beta*xb_1d(lc) +
     $          (1.-oprim(lc)+oprim(lc)*beta+zb_a(lc)*cmu)*xb_1u(lc))/
     $          denomb
           yb_1u(lc) = (  oprim(lc)*beta*xb_1u(lc) +
     $          (1.-oprim(lc)+oprim(lc)*beta-zb_a(lc)*cmu)*xb_1d(lc))/
     $          denomb
           z0p = xb_0u(lc) - cmu*yb_1d(lc)
           z0m = xb_0d(lc) + cmu*yb_1u(lc)
           yb_0d(lc) = (  oprim(lc)*beta*z0m +
     $          (1.-oprim(lc)+oprim(lc)*beta+zb_a(lc)*cmu)*z0p)/
     $          denomb
           yb_0u(lc) = (  oprim(lc)*beta*z0p +
     $          (1.-oprim(lc)+oprim(lc)*beta-zb_a(lc)*cmu)*z0m)/
     $          denomb
*     
         ENDIF
*
         IF ( planck  ) THEN
*
* Set coefficients in KST(22) for thermal source
*
* Calculate alpha coefficient 
*
            small = 1.E+20*smallest
            q0 = (1.-oprim(lc)) * pkag(lc-1)
            q1 = (1.-oprim(lc)) * pkagc(lc)
            q2 = (1.-oprim(lc)) * pkag(lc)
            deltat = taucpr(lc) - taucpr(lc-1)
*
* Case 1: source small at bottom layer
*
            IF ( (q2.LT.(q0*1.E-02) .OR. q2.LE.small )
     $           .AND. q1.GT.small .AND. q0.GT.small ) THEN
*
* alpha Eq. KS(50)
*
               zp_a(lc) = (2./deltat) * LOG( q0/q1 )
               IF ( ABS(zp_a(lc)) .GT. big )    zp_a(lc) = big
               IF ( ABS(zp_a(lc)*taucpr(lc-1)) .GE. ALOG(big) 
     $              .OR. ABS(zp_a(lc)*taucpr(lc)) .GE. ALOG(big) )
     $              THEN
                  xp_0(lc) =  big
               ELSE
                  xp_0(lc) = q0
               ENDIF
               xp_1(lc) = 0.0
*     
* Case 2: Source small at center and bottom of layer
*
            ELSE IF ( (q2.LE.(q1*1.E-02) .OR. q2.LE.small ) .AND.
     $              ((q1.LE.(q0*1.E-02)) .OR. (q1.LE.small))
     $              .AND. (q0.GT.small) ) THEN
*     
               zp_a(lc)  =   big / taucpr(ncut)
               xp_0(lc) = q0
               xp_1(lc) = 0.0
*     
*     Case 3:All sources zero
*
            ELSE IF ( q2.LE.small .AND. q1.LE.small
     $              .AND. q0.LE.small) THEN
               zp_a(lc)  = 0.0
               xp_0(lc) = 0.0
               xp_1(lc) = 0.0
*     
*     Case 4: Sources same at center, bottom and top of layer
*     or layer optically very thin
*
            ELSE IF ( (ABS((q2-q0)/q2).LT.1.E-04) .AND.
     $              (ABS((q2-q1)/q2).LT.1.E-04)
     $              .OR. deltat.LT. 1.E-04           ) THEN
*     
               zp_a(lc)  = 0.0
               xp_0(lc) = q0
               xp_1(lc) = 0.0
*     **  Case 5: Normal case
            ELSE
               arg = (q1/q2)**2. - q0/q2
               IF ( arg .LT. 0.0 ) arg = 0.0
*     
* alpha Eq. (44). For source that has its maximum at the top of the
* layer, use negative solution
*
               sgn = 1.0
               IF ( pkag(lc-1) .GT. pkag(lc) ) sgn = -1.
               fact3 = LOG(q1/q2 + sgn*SQRT(arg) )
               IF ( ABS(fact3) .LE. 0.005 ) THEN ! Be careful with log of
                  q1 = 0.99 * q1 ! numbers close to one
                  fact3 = LOG(q1/q2 + sgn*SQRT(arg) )
               ENDIF 
               zp_a(lc) = (2./deltat) * fact3
*aky               IF( ABS(zp_a(lc)*taucpr(lc)) .GT. 
*aky     $              (LOG(biggest) ) .OR. 
*aky     $              ABS(zp_a(lc)*taucpr(lc-1)) .GT. 
*aky     $              (LOG(biggest) )) THEN
               IF( ABS(zp_a(lc)*taucpr(lc)) .GT.  large .OR. 
     $              ABS(zp_a(lc)*taucpr(lc-1)) .GT. large ) THEN
                  zp_a(lc) = 0.0 
                  ENDIF
*     
* Dither alpha if it is close to an eigenvalue
*
               denomp =  fact1 * fact2 - (zp_a(lc)*cmu)*(zp_a(lc)*cmu)
               IF ( denomp .LT. 1.E-03 ) THEN
                  zp_a(lc) = 1.01*zp_a(lc)
               ENDIF
*
* Set constants in Eqs. KST(22)
*
               IF ( deltat .LT. 1.E-07 ) THEN
                  xp_1(lc) = 0.0
               ELSE
                  xp_1(lc) = (1./deltat)*(q2*EXP(zp_a(lc)*taucpr(lc)) 
     $                 -q0*EXP(zp_a(lc)*taucpr(lc-1)))
               ENDIF
               xp_0(lc) = q0 * EXP(zp_a(lc)*taucpr(lc-1)) -
     $              xp_1(lc)*taucpr(lc-1)
            ENDIF
*     
* Calculate particular solutions Eqs. KST(24-25) for internal thermal source
*
            denomp =  fact1 * fact2 - (zp_a(lc)*cmu)*(zp_a(lc)*cmu)
            yp_1d(lc) = (  oprim(lc)*beta*xp_1(lc) +
     $           (1.-oprim(lc)+oprim(lc)*beta+zp_a(lc)*cmu)*xp_1(lc))/
     $           denomp
            yp_1u(lc) = (  oprim(lc)*beta*xp_1(lc) +
     $           (1.-oprim(lc)+oprim(lc)*beta-zp_a(lc)*cmu)*xp_1(lc))/
     $           denomp
            z0p = xp_0(lc) - cmu*yp_1d(lc)
            z0m = xp_0(lc) + cmu*yp_1u(lc)
            yp_0d(lc) = (  oprim(lc)*beta*z0m +
     $           (1.-oprim(lc)+oprim(lc)*beta+zp_a(lc)*cmu)*z0p)/
     $           denomp
            yp_0u(lc) = (  oprim(lc)*beta*z0p +
     $           (1.-oprim(lc)+oprim(lc)*beta-zp_a(lc)*cmu)*z0m)/
     $           denomp
*     
         END IF
*     
 100  CONTINUE
*
* ===================  End loop on computational layers  ===============
*
      RETURN
      END
*
      SUBROUTINE settwo( biggest, bplank, btemp, ch, chtau, cmu,
     $     deltam, dtauc, dtaucpr, expbea, fbeam, flyr, gg, ggprim,
     $     layru, lyrcut, ncut, newgeo, nlyr, nn, nstr, ntau,
     $     oprim, pi, pkag, pkagc, planck, radius, smallest, spher, 
     $     ssalb, tauc, taucpr, temis, temper, tplank, ttemp, umu0,
     $     usrtau, utau, utaupr, wvnmlo, wvnmhi, zd)
*
* Perform miscellaneous setting-up operations
*
* Routines called:  errmsg, zeroit
*
* Input :  All are 'twostr' input variables (see doc file)
*
* Output:  ntau,utau   If usrtau = false
*          bplank      Intensity emitted from bottom boundary
*          ch          The Chapman factor
*          cmu         Computational polar angle
*          expbea      Transmission of direct beam
*          flyr        Truncated fraction in delta-m method
*          layru       Computational layer in which -utau- falls
*          lyrcut      Flag as to whether radiation will be zeroed
*                      below layer -ncut-
*          ncut        Computational layer where absorption
*                      optical depth first exceeds -abscut-
*          nn          nstr / 2  =  1
*          nstr        No. of streams (=2)
*          oprim       Delta-m-scaled single-scatter albedo
*          pkag,c      Planck function in each layer
*          taucpr      Delta-m-scaled optical depth
*          tplank      Intensity emitted from top boundary
*          utaupr      Delta-m-scaled version of -utau-
*
* Internal variables
*          abscut      Absorption optical depth, medium is cut off below
*                      this depth
*          tempc       Temperature at center of layer, assumed
*                      to be average of layer boundary temperatures
*
      DOUBLE PRECISION d1mach

c +-------------------------------------------------------------------+
c
c  LOCAL SYMBOLIC DIMENSIONS 
c
c  MXCLY  = Max no. of computational layers
c  MXULV  = Max no. of output levels
c  MXCMU  = Max no. of computation polar angles
c  MXUMU  = Max no. of output polar angles
c  MXPHI  = Max no. of output azimuthal angles
c +-------------------------------------------------------------------+

      INTEGER MXCLY, MXULV, MXTCMU, MI9M2, NNLYRI
      INCLUDE "DISORT.MXD"
      PARAMETER ( MXTCMU = 2, MI9M2 = 4,
     &          NNLYRI = MXTCMU*MXCLY)

      LOGICAL  deltam, lyrcut, newgeo, planck, spher, usrtau, first
     $     , hit
      INTEGER  layru(*)
      REAL ch(*), chtau(0:*), cmu, dtauc(*), dtaucpr(*),
     $     expbea(0:*), flyr(*), gg(*), ggprim(*), oprim(*), pkag(0:*),
     $     pkagc(*), ssalb(*), tauc(0:*), taucpr(0:*),
     $     temper(0:*), utau(*), utaupr(*), zd(0:*)
      INTEGER nfac1(mxcly), nfac2(mxcly)
      REAL fac1(mxcly,2*mxcly), fac2(mxcly,2*mxcly)
      REAL chapman
*
      SAVE nfac1, nfac2
      SAVE first, abscut
      DATA  abscut / 10. /, first / .TRUE. /
*
      IF ( first ) THEN
         first    = .FALSE.
         smallest = r1mach(1)
         biggest  = r1mach(2)
         nstr     = 2
         nn       = nstr / 2
      ENDIF
*
      IF ( .NOT.usrtau ) THEN
*
* Set output levels at computational layer boundaries
*
         ntau = nlyr + 1
         DO 30  lc = 0, ntau-1
            utau(lc+1) = tauc(lc)
 30      CONTINUE
      END IF
*
* Apply delta-m scaling and move description of computational
* layers to local variables
*
      CALL  zeroit( taucpr(0), mxcly+1 )
      CALL  zeroit( expbea(0), mxcly+1 )
      CALL  zeroit( flyr, mxcly )
      CALL  zeroit( oprim, mxcly )
      abstau = 0.0
      DO  60  lc = 1, nlyr
         IF ( abstau.LT.abscut )  ncut = lc
         abstau = abstau + ( 1. - ssalb(lc) ) * dtauc(lc)
*
         IF ( .NOT.deltam )  THEN
            oprim(lc)  = ssalb(lc)
            taucpr(lc) = tauc(lc)
            f          = 0.0
            ggprim(lc) = gg(lc)
            dtaucpr(lc)= dtauc(lc)
         ELSE
*
* Do delta-m transformation Eqs. WW(20a,20b,14)
*
            f = gg(lc) * gg(lc)
            taucpr(lc) = taucpr(lc-1) + ( 1. - f*ssalb(lc) ) * dtauc(lc)
            oprim(lc)  = ssalb(lc) * ( 1. - f ) / ( 1. - f * ssalb(lc) )
            ggprim(lc) =  (gg(lc)-f) / (1.-f)
            dtaucpr(lc)= taucpr(lc) - taucpr(lc-1)
         ENDIF
*
         flyr(lc)   = f
*
 60   CONTINUE
*
* If no thermal emission, cut off medium below absorption optical
* depth = abscut ( note that delta-m transformation leaves absorption
* optical depth invariant ).  Not worth the trouble for one-layer
* problems, though.
*
      lyrcut = .FALSE.
      IF ( abstau.GE.abscut .AND. .NOT. planck
     $     .AND. nlyr.GT.1 )  lyrcut =.TRUE.
      IF ( .NOT.lyrcut )  ncut = nlyr
*
* Calculate chapman function is spherical geometry, set expbea and ch
* for beam source.
*
      IF ( fbeam.GT.0.0 )  THEN
         zenang      = ACOS(umu0) * 180. / pi
         IF ( spher .AND. newgeo ) THEN
            z_lay = 0.0
            CALL geofac( fac1, nfac1, z_lay, nlyr, zd, zenang, radius )
            z_lay = 0.5
            CALL geofac( fac2, nfac2, z_lay, nlyr, zd, zenang, radius )
         ENDIF
         chtau(0) = 0.0
         IF ( spher ) THEN
            expbea( 0 ) = 1.0
            IF( umu0 .LT. 0.0 ) THEN
               expbea(0) =
     $              EXP(-chapman( 1, nfac1, fac1, dtaucpr ))
            ENDIF
            DO lc = 1, ncut
               chtau(lc) = chapman( lc, nfac1, fac1, dtaucpr )
               taup       = tauc(lc-1) + dtauc(lc)/2.0
               ch(lc)     = taup/chapman( lc, nfac2, fac2, dtauc )
               expbea(lc) = EXP(-chapman( lc, nfac1, fac1, dtaucpr ) )
            ENDDO
         ELSE IF ( .NOT. spher ) THEN
            expbea( 0 ) = 1.0
            DO lc = 1, ncut
               ch(lc)     = umu0
               chtau(lc)  = taucpr(lc) / umu0 
               expbea(lc) = EXP( - taucpr(lc) / umu0 )
            ENDDO
         ENDIF
      ENDIF
*
* Set arrays defining location of user output levels within 
* delta-m-scaled computational mesh
* 
      DO  lu = 1, ntau
         hit = .FALSE.
         DO  lc = 1, nlyr
            IF ( utau(lu).GE.tauc(lc-1) .AND. utau(lu).LE.tauc(lc) 
     $           .AND. .NOT. hit) THEN
               utaupr(lu) = utau(lu)
               IF(deltam) utaupr(lu) = taucpr(lc-1) + (1.-ssalb(lc)
     $              * flyr(lc)) * (utau(lu) - tauc(lc-1))
               layru(lu) = lc
               hit = .TRUE.
            ENDIF
         ENDDO
         IF ( .NOT. hit ) THEN
            lc = nlyr
            utaupr(lu) = utau(lu)
            IF(deltam) utaupr(lu) = taucpr(lc-1) + (1.-ssalb(lc)
     $           * flyr(lc)) * (utau(lu) - tauc(lc-1))
            layru(lu) = lc
            hit = .TRUE.
         ENDIF
      ENDDO
*
* Set computational polar angle cosine for double gaussian 
* quadrature; cmu = 0.5, or  single gaussian quadrature; cmu = 1./sqrt(3)
* See KST for discussion of which is better for your specific application
*
      IF ( planck .AND. fbeam .EQ. 0.0 ) THEN
         cmu =  0.5
      ELSE
         cmu = 1./SQRT(3.0)
      ENDIF
*
* Calculate planck functions
*
      IF ( .NOT. planck )  THEN
         bplank = 0.0
         tplank = 0.0
         CALL  zeroit( pkag, mxcly+1 )
         CALL  zeroit( pkagc, mxcly )
      ELSE
         tplank = temis * tplkavg( wvnmlo, wvnmhi, ttemp )
         bplank =         tplkavg( wvnmlo, wvnmhi, btemp )
         DO 180  lev = 0, nlyr
            pkag( lev ) = tplkavg( wvnmlo, wvnmhi, temper(lev) )
 180     CONTINUE
         DO 190 lc = 1, nlyr
            tempc = 0.5 * ( temper(lc-1) + temper(lc) )
            pkagc( lc ) = tplkavg( wvnmlo, wvnmhi, tempc )
 190     CONTINUE
      END IF
      RETURN
      END
*
      SUBROUTINE solbcd( albedo, b, bplank, cband, cmu, diag, expbea,
     $     fbeam, fisot, ipvt, kk, ll, lyrcut, mi9m2, mxtcmu,
     $     ncut, nnlyri, pi, rr, subd,
     $     superd, sol, tplank, taucpr, umu0, yb_0d, yb_0u, yb_1d,
     $     yb_1u, yp_0d, yp_0u, yp_1d, yp_1u, zb_a, zp_a )
*
* Construct right-hand side vector -b- for general boundary conditions
* and solve system of equations obtained from the boundary conditions
* and the continuity-of-intensity-at-layer-interface equations.
*
* Routines called: sgbfa, sgbsl, zeroit
*
* I n p u t      v a r i a b l e s:
*
*       bplank   :  Bottom boundary thermal emission
*       cband    :  Left-hand side matrix of linear system Eqs. KST(38-41)
*                   in banded form required by linpack solution routines
*       cmu      :  Abscissae for gauss quadrature over angle cosine
*       expbea   :  Transmission of incident beam, EXP(-taucpr/ch)
*       lyrcut   :  Logical flag for truncation of comput. layer
*       ncut     :  Total number of computational layers considered
*       tplank   :  Top boundary thermal emission
*       taucpr   :  Cumulative optical depth (delta-m-scaled)
*       yb_0d,u, :  Thermal source vectors, Eq. KST(24-25)
*        yb_1d,u     
*       yp_0d,u, :  Beam source vectors, Eq. KST(24-25)
*        yp_1d,u     
*       zb_a     :  Beam source coefficient alfa in Eq. KST(22)
*       zp_a     :  Thermal source coefficient alfa in Eq. KST(22)
*       (Remainder are 'twostr' input variables)
*
* O u t p u t     v a r i a b l e s:
*
*       b        :  Right-hand side vector of Eqs. KST(38-41) going into
*                   *sgbsl*; returns as solution vector of Eqs. KST(38-41)
*                   constants of integration 
*      ll        :  Permanent storage for -b-, but re-ordered
*
* I n t e r n a l    v a r i a b l e s:
*
*       it       :  Pointer for position in  -b-
*       ncd      :  Number of diagonals below or above main diagonal
*       rcond    :  Indicator of singularity for -cband-
*       z        :  Scratch array required by *sgbco*
*+---------------------------------------------------------------------+
*
      INTEGER  ipvt(*), itrierr
      LOGICAL  lyrcut
      REAL b(*), cband( mi9m2,nnlyri ), cmu, diag(*), expbea(0:*), 
     $     kk(*), ll( mxtcmu,* ), rr(*), subd(*), superd(*), sol(*),
     $     taucpr( 0:* ), yb_0d(*), yb_0u(*), yb_1d(*), yb_1u(*),
     $     yp_0d(*), yp_0u(*), yp_1d(*), yp_1u(*), zb_a(*), zp_a(*)
      REAL gam2(252)
* 
* First top row, top boundary condition
*
      irow        = 1
      lc          = 1
*     subd(irow)  = undefined
      diag(irow)  = rr(lc) * EXP(-kk(lc) * taucpr(lc))
      superd(irow)= 1.0
*aky      diag(irow)  = 1.0
*aky      superd(irow)= rr(lc) * EXP(-kk(lc) * taucpr(lc))
*aky      write(*,*)"diag", irow, diag(irow), rr(lc),kk(lc),taucpr(lc)
*
* next from layer no. 2 to nlyr -1
*
      nloop = ncut - 1 
      DO lc = 1, nloop
         irow         = irow + 1
         wk0          = EXP(-kk(lc) * (taucpr(lc) - taucpr(lc-1)))
         wk1          = EXP(-kk(lc+1) * (taucpr(lc+1) - taucpr(lc)))
         subd(irow)   = 1.0 - rr(lc) * rr(lc+1)
         diag(irow)   = ( rr(lc) -  rr(lc+1 ) ) * wk0
*aky         subd(irow)   = ( rr(lc) -  rr(lc+1 ) ) * wk0
*aky         diag(irow)   = 1.0 - rr(lc) * rr(lc+1)
         superd(irow) = - ( 1. - rr(lc+1) * rr(lc+1 ) ) * wk1
         irow         = irow + 1
         subd(irow)   = ( 1.0 - rr(lc) * rr(lc) ) * wk0 
         diag(irow)   = ( rr(lc) -  rr(lc+1 ) ) * wk1
*aky         subd(irow)   = ( rr(lc) -  rr(lc+1 ) ) * wk1
*aky         diag(irow)   = ( 1.0 - rr(lc) * rr(lc) ) * wk0 
         superd(irow) = - ( 1. - rr(lc+1) * rr(lc ) ) 
*aky      write(*,*)"diag", irow, subd(irow), diag(irow), superd(irow)
      ENDDO
*
* bottom layer
*
      irow         = irow + 1
      lc           = ncut
*     superd(irow) = undefined
      wk           = EXP( -kk(lc) * (taucpr(lc) - taucpr(lc-1)) )
      IF ( lyrcut ) THEN
         subd(irow) = 1.0
         diag(irow) = rr(lc) * wk
*aky         subd(irow) = rr(lc) * wk
*aky         diag(irow) = 1.0
      ELSE
         subd(irow) = 1.0 - 2.0*albedo*cmu*rr(lc)
         diag(irow) = ( rr(lc) - 2.0*albedo*cmu ) * wk         
*aky         subd(irow) = ( rr(lc) - 2.0*albedo*cmu ) * wk         
*aky         diag(irow) = 1.0 - 2.0*albedo*cmu*rr(lc)
      ENDIF
*
      CALL  zeroit( b, nnlyri )
*
* Construct -b-,  for parallel beam + bottom reflection +
* thermal emission at top and/or bottom
*
* Top boundary, right-hand-side of Eq. KST(28)
*
      lc      = 1
      irow    = 1
      b(irow) = - yb_0d(lc) - yp_0d(lc) +fisot +tplank
*
* Continuity condition for layer interfaces,
* right-hand-side of Eq. KST(29)
*
      DO   lc = 1, nloop
         rpp1_m = EXP(-zb_a(lc+1)*taucpr(lc))*
     $        (yb_0d(lc+1)+yb_1d(lc+1)*taucpr(lc)) 
     $        + EXP(-zp_a(lc+1)*taucpr(lc))*
     $        (yp_0d(lc+1)+yp_1d(lc+1)*taucpr(lc))
*aky         write(*,*)"rpp1_m", lc, zp_a(lc+1)*taucpr(lc),
*aky     $        yp_0d(lc+1),yp_1d(lc+1), EXP(-zp_a(lc+1)*taucpr(lc))
         rp_m   =  EXP(-zb_a(lc)*taucpr(lc))*
     $        (yb_0d(lc)+yb_1d(lc)*taucpr(lc))
     $        +  EXP(-zp_a(lc)*taucpr(lc))*
     $        (yp_0d(lc)+yp_1d(lc)*taucpr(lc))
         rpp1_p = EXP(-zb_a(lc+1)*taucpr(lc))*
     $        (yb_0u(lc+1)+yb_1u(lc+1)*taucpr(lc))
     $        +  EXP(-zp_a(lc+1)*taucpr(lc))*
     $        (yp_0u(lc+1)+yp_1u(lc+1)*taucpr(lc))
         rp_p   = EXP(-zb_a(lc)*taucpr(lc))*
     $        (yb_0u(lc)+yb_1u(lc)*taucpr(lc))
     $        +  EXP(-zp_a(lc)*taucpr(lc))*
     $        (yp_0u(lc)+yp_1u(lc)*taucpr(lc))
         irow    = irow + 1
         b(irow) = rpp1_p - rp_p - rr(lc+1) * ( rpp1_m - rp_m )
         irow    = irow + 1
         b(irow) = rpp1_m - rp_m - rr(lc) * ( rpp1_p - rp_p )
*aky         write(*,*)"bb", irow, rpp1_p, rp_p, rpp1_m, rp_m
      ENDDO
*
* Bottom boundary
*
      irow = irow + 1
      lc   = ncut
      IF ( lyrcut ) THEN
*
* Right-hand-side of Eq. KST(30)
*
         b(irow) = - EXP(-zb_a(ncut)*taucpr(ncut))*
     $        (yb_0u(ncut)+yb_1u(ncut)*taucpr(ncut))
     $        - EXP(-zp_a(ncut)*taucpr(ncut))*
     $        (yp_0u(ncut)+yp_1u(ncut)*taucpr(ncut))
      ELSE
         sum = cmu * albedo*(EXP(-zb_a(ncut)*taucpr(ncut))*
     $        (yb_0d(ncut)+yb_1d(ncut)*taucpr(ncut))
     $        +  EXP(-zp_a(ncut)*taucpr(ncut))*
     $        (yp_0d(ncut)+yp_1d(ncut)*taucpr(ncut)))
         refflx = 1.
         IF ( umu0 .LE. 0.0 ) refflx = 0.0
         b(irow) = 2.*sum + 
     $        ( albedo * umu0*fbeam/pi*refflx ) * expbea(ncut) 
     $        + (1.-albedo) * bplank
     $        -  EXP(-zb_a(ncut)*taucpr(ncut))*
     $        (yb_0u(ncut)+yb_1u(ncut)*taucpr(ncut))
     $        -  EXP(-zp_a(ncut)*taucpr(ncut))*
     $        (yp_0u(ncut)+yp_1u(ncut)*taucpr(ncut))
      END IF
*
* solve for constants of integration by inverting matrix KST(38-41)
*
      nrow = irow
*
* If you want to use a simple tridiagonal solver instead of the 
* provided linpack routines, comment the lines contained within
* * *linpack stuff* and uncomment the next(+1) eight lines.
*
*aky      write(*,*) "before tridag", nrow
*aky      CALL tridag( subd, diag, superd, b, sol, gam2, nrow, itrierr )
*aky      write(*,*) "finished tridag", itrierr
*aky      irow = 1
*aky      DO lc = 1, ncut
*aky         ll(1,lc) = sol(irow)      ! downward direction
*aky         write(*,*) "ll", lc, irow, irow+1, sol(irow), sol(irow+1)
*aky         irow     = irow + 1
*aky         ll(2,lc) = sol(irow)      ! upward direction
*aky         irow     = irow + 1
*aky      ENDDO
*
* *linpack stuff*
*
          CALL zeroit( cband, mi9m2*nnlyri)
          DO irow = 1, nrow
             cband(1,irow)   = 0.0
             cband(3,irow)   = diag(irow)
          ENDDO
          DO irow = 1, nrow-1
             cband(2,irow+1) = superd(irow)
          ENDDO
          DO irow = 2, nrow
             cband(4,irow-1) = subd(irow)
          ENDDO
*
          CALL sgbfa(cband, mi9m2, nrow, 1, 1, ipvt, info )
          job = 0
          CALL sgbsl(cband, mi9m2, nrow, 1, 1, ipvt, b, job )
*
* unpack
*
          irow = 1
          DO lc = 1, ncut
             ll(1,lc) = b(irow)      ! downward direction
             irow     = irow + 1
             ll(2,lc) = b(irow)      ! upward direction
             irow     = irow + 1
          ENDDO
*
* *linpack stuff*
*
      RETURN
      END
*
      SUBROUTINE  tchekin(albedo, btemp, dtauc, fbeam, fisot, gg,
     $     ierror, maxcly, maxulv, mxcly, mxulv, nlyr, ntau, planck,
     $     quiet, spher, ssalb, tauc, temis, temper, ttemp, umu0,
     $     usrtau, utau, wvnmlo, wvnmhi, zd )
*
* Checks the input dimensions and variables
*
      LOGICAL inperr, planck, quiet, spher, usrtau, wrtbad, wrtdim
      LOGICAL TEMPSTEPWARN
      INTEGER  maxcly, maxulv, mxcly, mxulv, nlyr, ntau, ierror(*)
      REAL albedo, btemp, dtauc(*), fbeam, fisot, gg(*), 
     $     ssalb(*), tauc(0:*), temis, temper(0:*), ttemp,
     $     umu0, utau(*), wvnmlo, wvnmhi, zd(0:*)
*
      DATA      TEMPSTEPWARN / .True. /
      SAVE      TEMPSTEPWARN

      inperr = .FALSE.
      IF ( nlyr.LT.1 ) THEN
         inperr = wrtbad( 'nlyr' )
         ierror(1) = 1
      ENDIF
      IF ( nlyr.GT.maxcly ) THEN
         inperr = wrtbad( 'maxcly' )
         ierror(2) = 1
      ENDIF
*
      DO 10  lc = 1, nlyr
         IF ( dtauc(lc).LT.0.0 )  THEN
            inperr = wrtbad( 'dtauc' )
            ierror(3) = ierror(3) + 1
         ENDIF
         IF ( ssalb(lc).LT.0.0 .OR. ssalb(lc).GT.1.0 ) THEN
            inperr = wrtbad( 'ssalb' )
            ierror(4) = ierror(4) + 1
         ENDIF
         IF ( planck )  THEN
            IF( lc.EQ.1 .AND. temper(0).LT.0.0 ) THEN
               inperr = wrtbad( 'temper' )
               ierror(5) = ierror(5) + 1
            ENDIF
            IF( temper(lc).LT.0.0 ) THEN 
               inperr = wrtbad( 'temper' )
               ierror(5) = ierror(5) + 1
            ENDIF
         ENDIF
         IF( gg(lc).LT.-1.0 .OR. gg(lc).GT.1.0 ) THEN
            inperr = wrtbad( 'gg' )
            ierror(6) = ierror(6) + 1
         ENDIF
 10   CONTINUE
      IF ( spher ) THEN
         DO 11 lc = 1, nlyr
            IF ( zd(lc) .GT. zd(lc-1) ) THEN
               inperr = wrtbad( 'zd' )
               ierror(7) = ierror(7) + 1
            ENDIF
 11      CONTINUE   
      ENDIF
*
      IF ( usrtau )  THEN
         IF ( ntau.LT.1 ) THEN
            inperr = wrtbad( 'ntau' )
            ierror(8) = 1
         ENDIF
         IF ( maxulv.LT.ntau ) THEN
            inperr = wrtbad( 'maxulv' )
            ierror(9) = 1
         ENDIF
         DO 20  lu = 1, ntau
            IF( ABS(utau(lu)-tauc(nlyr) ).LE.1.E-4 .OR.
     $           (ABS(utau(lu)-tauc(nlyr)).LE.1.E-6*tauc(nlyr)))
     $           utau(lu) = tauc(nlyr)
            IF( utau(lu).LT.0.0 .OR. utau(lu).GT.tauc(nlyr) ) THEN
               inperr = wrtbad( 'utau' )
               ierror(10) = ierror(10) + 1
            ENDIF
 20      CONTINUE
      ELSE
         IF ( maxulv.LT.nlyr+1 ) THEN
            inperr = wrtbad( 'maxulv' )
            ierror(11) = 1
         ENDIF
      END IF
*
      IF ( fbeam.LT.0.0 ) THEN
         inperr = wrtbad( 'fbeam' )
         ierror(12) = 1
      ENDIF
*
      IF ( .NOT. spher ) THEN
         umumin = 0.0
         IF ( fbeam.GT.0.0 .AND. ( umu0.LT. umumin .OR. umu0.GT.1.0 ) )
     $        THEN
            inperr = wrtbad( 'umu0' )
            ierror(13) = 1
         ENDIF
      ELSE IF ( spher ) THEN
         umumin = - 1.0
         IF ( fbeam.GT.0.0 .AND. ( umu0.LE.umumin .OR. umu0.GT.1.0 ) )
     $        THEN
            inperr = wrtbad( 'umu0' )
            ierror(13) = 1
         ENDIF
      ENDIF
*
      IF ( fisot.LT.0.0 ) THEN
         inperr = wrtbad( 'fisot' )
         ierror(14) = 1
      ENDIF
*
      IF ( albedo.LT.0.0 .OR. albedo.GT.1.0 ) THEN
         inperr = wrtbad( 'albedo' )
         ierror(15) = 1
      ENDIF
*
      IF ( planck )  THEN
         IF ( wvnmlo.LT.0.0 .OR. wvnmhi.LT.wvnmlo ) THEN
            inperr = wrtbad( 'wvnmlo,hi' )
            ierror(16) = 1
         ENDIF
         IF ( temis.LT.0.0 .OR. temis.GT.1.0 ) THEN
            inperr = wrtbad( 'temis' )
            ierror(17) = 1
         ENDIF
         IF ( btemp.LT.0.0 ) THEN
            inperr = wrtbad( 'btemp' )
            ierror(18) = 1
         ENDIF
         IF ( ttemp.LT.0.0 ) THEN
            inperr = wrtbad( 'ttemp' )
            ierror(19) = 1
         ENDIF
      END IF
*
      IF ( mxcly.LT.nlyr ) THEN
         inperr = wrtdim( 'mxcly', nlyr )
         ierror(20) = 1
      ENDIF
      IF ( usrtau .AND. mxulv.LT.ntau ) THEN
         inperr = wrtdim( 'mxulv', ntau )
         ierror(21) = 1
      ENDIF
      IF ( .NOT.usrtau .AND. mxulv.LT.nlyr+1 ) THEN
         inperr = wrtdim( 'mxulv', nlyr+1 )
         ierror(22) = 1
      ENDIF
*
      IF ( inperr )
     $   CALL errmsg( 'twostr--input and/or dimension errors', .TRUE. )
*
      DO 100  lc = 1, nlyr
         IF ( (planck .AND. ABS(temper(lc)-temper(lc-1)) .GT. 50.0) 
     $          .AND. .NOT. quiet .AND. TEMPSTEPWARN) THEN
            CALL errmsg( 'chekin--vertical temperature step larger'
     $           // ' than 50K may be too large for good accuracy',
     $               .FALSE. )
            WRITE (0,*) '   T[',lc-1,']=',temper(lc-1),
     $                   ', T[',lc  ,']=',temper(lc)
            WRITE (0,*) ''
            TEMPSTEPWARN = .FALSE.
         ENDIF

100   CONTINUE
*
      RETURN
      END
*
      SUBROUTINE tfluxes( ch, cmu, dfdt, fbeam, flup, fldn, fldir, kk,
     $     layru, ll, lyrcut, maxulv, mxtcmu, mxulv, ncut, ntau,
     $     pi, planck, prnt, oprim, rfldir, rfldn, rr, spher, ssalb,
     $     taucpr, u0c, uavg, umu0, utau, utaupr, xp_0, xp_1,
     $     yb_0d, yb_0u, yb_1d, yb_1u, yp_0d, yp_0u, yp_1d, yp_1u,
     $     zb_a, zp_a )
*
* Calculates the radiative fluxes, mean intensity, and flux
* derivative with respect to optical depth from the 
* azimuthally-averaged intensity
*
* I n p u t     v a r i a b l e s:
*
*       ch       :  Chapman factor
*       cmu      :  Abscissae for gauss quadrature over angle cosine
*       kk       :  Eigenvalues 
*       layru    :  Layer number of user level -utau-
*       ll       :  Constants of integration in Eqs. KST(42-43), obtained
*                   by solving Eqs. KST(38-41)
*       lyrcut   :  Logical flag for truncation of comput. layer
*       ncut     :  Number of computational layer where absorption
*                     optical depth exceeds -abscut-
*       oprim    :  Delta-m scaled single scattering albedo
*       rr       :  Eigenvectors at polar quadrature angles 
*       taucpr   :  Cumulative optical depth (delta-m-scaled)
*       utaupr   :  Optical depths of user output levels in delta-m
*                     coordinates;  equal to  -utau- if no delta-m
*       xp_0,    :  Thermal source coefficients x-sub-zero and
*        xp_1         x-sub-one in Eq. KST(22)
*       yb_0d,u, :  Thermal source vectors, Eq. KST(23)
*        yb_1d,u     
*       yp_0d,u, :  Beam source vectors, Eq. KST(23)
*        yp_1d,u     
*       zb_a     :  Beam source coefficient alfa in Eq. KST(22)
*       zp_a     :  Thermal source coefficient alfa in Eq. KST(22)
*       (remainder are 'twostr' input variables)
*
* O u t p u t     v a r i a b l e s:
*
*       u0c      :  Azimuthally averaged intensities at polar
*                     quadrature angle cmu
*       (rfldir, rfldn, flup, dfdt, uavg are 'twostr' output variables)
*
* I n t e r n a l       v a r i a b l e s:
*
*       dirint   :  direct intensity attenuated
*       fdntot   :  total downward flux (direct + diffuse)
*       fldir    :  direct-beam flux (delta-m scaled)
*       fldn     :  diffuse down-flux (delta-m scaled)
*       fnet     :  net flux (total-down - diffuse-up)
*       fact     :  EXP( - utaupr / ch ), where ch is the Chapman factor
*       plsorc   :  Planck source function (thermal)
*+---------------------------------------------------------------------+
      INTEGER layru(*)
      LOGICAL lyrcut, planck, prnt(*), spher
      REAL ch(*), cmu, dfdt(*), flup(*), fldir(*), fldn(*), kk( * ), 
     $     ll( mxtcmu,* ), oprim(*), rfldir(*), rfldn(* ), rr(*), 
     $     ssalb(*), taucpr( 0:* ), u0c( mxtcmu,mxulv ), uavg(*),
     $     utau(*), utaupr(*), xp_0(*), xp_1(*), yb_0d(*), yb_0u(*),
     $     yb_1d(*), yb_1u(*), yp_0d(*), yp_0u(*), yp_1d(*), yp_1u(*),
     $     zb_a(*), zp_a(*)
*
      IF ( prnt(2) )  WRITE( *,1010 )
*
* Zero twostr output arrays
*
      CALL  zeroit( u0c, mxulv*mxtcmu )
      CALL  zeroit( rfldir, maxulv )
      CALL  zeroit( fldir,  mxulv )
      CALL  zeroit( rfldn,  maxulv )
      CALL  zeroit( fldn,   mxulv )
      CALL  zeroit( flup,   maxulv )
      CALL  zeroit( uavg,   maxulv )
      CALL  zeroit( dfdt,   maxulv )
*
* Loop over user levels
*
      IF ( planck ) THEN
         DO 100  lu = 1, ntau
            lyu = layru(lu)
            iq = 1
            u0c( iq,lu ) = u0c( iq,lu ) +
     $           EXP(-zp_a(lyu)*utaupr(lu))*
     $           (yp_0d(lyu)+yp_1d(lyu)*utaupr(lu))
            iq = 2
            u0c( iq,lu ) = u0c( iq,lu ) + 
     $           EXP(-zp_a(lyu)*utaupr(lu))*
     $           (yp_0u(lyu)+yp_1u(lyu)*utaupr(lu))
 100     CONTINUE
      ENDIF
*
      DO 102  lu = 1, ntau
*
         lyu = layru(lu)
         IF ( lyrcut .AND. lyu.GT.ncut ) THEN
*
* No radiation reaches this level
*
            fdntot = 0.0
            fnet   = 0.0
            plsorc = 0.0
            GO TO 90                    ! ** Done with this layer
         END IF
*
         IF ( fbeam.GT.0.0 ) THEN
            iq = 1
            u0c( iq,lu ) = u0c(iq,lu) + 
     $           EXP(-zb_a(lyu)*utaupr(lu)) *
     $           (yb_0d(lyu)+yb_1d(lyu)*utaupr(lu) )
            iq = 2
            u0c( iq,lu ) = u0c(iq,lu) +
     $           EXP(-zb_a(lyu)*utaupr(lu))*
     $           ( yb_0u(lyu)+yb_1u(lyu)*utaupr(lu) )
            IF ( umu0.GT.0.0 .OR. spher  ) THEN
               fact         = EXP( - utaupr(lu) / ch(lyu) )
               dirint       = fbeam * fact
               fldir(  lu ) = ABS(umu0) * ( fbeam * fact )
               rfldir( lu ) = ABS(umu0)*fbeam*EXP( -utau(lu)/ch(lyu) )

            ELSE
               dirint       = 0.0
               fldir(  lu ) = 0.0
               rfldir( lu ) = 0.0
            ENDIF
         ELSE
            dirint       = 0.0
            fldir(  lu ) = 0.0
            rfldir( lu ) = 0.0
         END IF
*
         iq = 1
         u0c( iq,lu ) =  u0c( iq,lu ) + rr(lyu) * ll(1,lyu) *
     $                EXP( - kk(lyu) * (taucpr(lyu) - utaupr(lu)) )
     $              + ll(2,lyu) *
     $                EXP( - kk(lyu) * (utaupr(lu) - taucpr(lyu-1)) )
         iq = 2
         u0c(iq,lu ) = u0c( iq,lu ) + ll(1,lyu) *
     $                EXP( - kk(lyu) * (taucpr(lyu) - utaupr(lu) ) )
     $              + rr(lyu) * ll(2,lyu) *
     $                EXP( - kk(lyu) * (utaupr(lu) - taucpr(lyu-1)) )
*
* Calculate fluxes and mean intensities
*
* Downward and upward fluxes from Eq. KST(9)
*
         fldn( lu )  = 2.0 * pi * cmu * u0c(1,lu)
         flup( lu )  = 2.0 * pi * cmu * u0c( 2,lu )
         fdntot      = fldn( lu ) + fldir( lu )
         fnet        = fdntot - flup( lu )
         rfldn( lu ) = fdntot - rfldir( lu )
*
* Mean intensity from Eq. KST(10)
*
         uavg(lu)   = u0c( 1,lu ) + u0c( 2, lu )
         uavg( lu ) = ( 2.0 * pi * uavg(lu) + dirint ) / ( 4.*pi )
*
* Flux divergence from Eqs. KST(11-12)
*
         plsorc = (1./(1.-oprim(lyu)))*EXP(-zp_a(lyu)*utaupr(lu))*
     $                          (xp_0(lyu) + xp_1(lyu)* utaupr(lu))
         dfdt( lu ) =  (1.0-ssalb(lyu)) * 4.*pi* (  uavg(lu) - plsorc )
*
 90      IF( prnt(2) )  WRITE( *,1020 ) utau(lu), lyu, rfldir(lu),
     $                                 rfldn(lu), fdntot, flup(lu),
     $                                 fnet, uavg(lu), plsorc, dfdt(lu)
 102  CONTINUE
*
 1010 FORMAT( //, 21X,
     $ '<----------------------- Fluxes ----------------------->', /,
     $ '   optical  compu    downward    downward    downward     ',
     $ ' upward                    mean      Planck   d(net flux)', /,
     $ '     depth  layer      direct     diffuse       total     ',
     $ 'diffuse         net   intensity      source   / d(op dep)', / )
 1020 FORMAT( F10.4, I7, 1P,7E12.3, E14.3 )
*
      RETURN
      END
*
      REAL FUNCTION  TPLKAVG ( WNUMLO, WNUMHI, T )
C
C        COMPUTES PLANCK FUNCTION INTEGRATED BETWEEN TWO WAVENUMBERS,
*        except if wnumlo .EQ. wnmuhi, then the Planck function at 
*        wnumlo is returned
C
C  NOTE ** CHANGE 'R1MACH' TO 'D1MACH' TO RUN IN DOUBLE PRECISION
C
C  I N P U T :  WNUMLO : LOWER WAVENUMBER ( INV CM ) OF SPECTRAL
C                           INTERVAL
C               WNUMHI : UPPER WAVENUMBER
C               T      : TEMPERATURE (K)
C
C  O U T P U T :  PLKAVG : INTEGRATED PLANCK FUNCTION ( WATTS/SQ M )
C                           = INTEGRAL (WNUMLO TO WNUMHI) OF
C                              2H C**2  NU**3 / ( EXP(HC NU/KT) - 1)
C                              (WHERE H=PLANCKS CONSTANT, C=SPEED OF
C                              LIGHT, NU=WAVENUMBER, T=TEMPERATURE,
C                              AND K = BOLTZMANN CONSTANT)
C
C  REFERENCE : SPECIFICATIONS OF THE PHYSICAL WORLD: NEW VALUE
C                 OF THE FUNDAMENTAL CONSTANTS, DIMENSIONS/N.B.S.,
C                 JAN. 1974
C
C  METHOD :  FOR  -WNUMLO-  CLOSE TO  -WNUMHI-, A SIMPSON-RULE
C            QUADRATURE IS DONE TO AVOID ILL-CONDITIONING; OTHERWISE
C
C            (1)  FOR WAVENUMBER (WNUMLO OR WNUMHI) SMALL,
C                 INTEGRAL(0 TO WNUM) IS CALCULATED BY EXPANDING
C                 THE INTEGRAND IN A POWER SERIES AND INTEGRATING
C                 TERM BY TERM;
C
C            (2)  OTHERWISE, INTEGRAL(WNUMLO/HI TO INFINITY) IS
C                 CALCULATED BY EXPANDING THE DENOMINATOR OF THE
C                 INTEGRAND IN POWERS OF THE EXPONENTIAL AND
C                 INTEGRATING TERM BY TERM.
C
C  ACCURACY :  AT LEAST 6 SIGNIFICANT DIGITS, ASSUMING THE
C              PHYSICAL CONSTANTS ARE INFINITELY ACCURATE
C
C  ERRORS WHICH ARE NOT TRAPPED:
C
C      * POWER OR EXPONENTIAL SERIES MAY UNDERFLOW, GIVING NO
C        SIGNIFICANT DIGITS.  THIS MAY OR MAY NOT BE OF CONCERN,
C        DEPENDING ON THE APPLICATION.
C
C      * SIMPSON-RULE SPECIAL CASE IS SKIPPED WHEN DENOMINATOR OF
C        INTEGRAND WILL CAUSE OVERFLOW.  IN THAT CASE THE NORMAL
C        PROCEDURE IS USED, WHICH MAY BE INACCURATE IF THE
C        WAVENUMBER LIMITS (WNUMLO, WNUMHI) ARE CLOSE TOGETHER.
C ----------------------------------------------------------------------
C                                   *** ARGUMENTS
      REAL     T, WNUMLO, WNUMHI
C                                   *** LOCAL VARIABLES
C
C        A1,2,... :  POWER SERIES COEFFICIENTS
C        c1       :  First radiation constant ( 2 * h * c**2 )
C        C2       :  H * C / K, IN UNITS CM*K (H = PLANCKS CONSTANT,
C                      C = SPEED OF LIGHT, K = BOLTZMANN CONSTANT)
C        D(I)     :  EXPONENTIAL SERIES EXPANSION OF INTEGRAL OF
C                       PLANCK FUNCTION FROM WNUMLO (I=1) OR WNUMHI
C                       (I=2) TO INFINITY
C        EPSIL    :  SMALLEST NUMBER SUCH THAT 1+EPSIL .GT. 1 ON
C                       COMPUTER
C        EX       :  EXP( - V(I) )
C        EXM      :  EX**M
C        MMAX     :  NO. OF TERMS TO TAKE IN EXPONENTIAL SERIES
C        MV       :  MULTIPLES OF 'V(I)'
C        P(I)     :  POWER SERIES EXPANSION OF INTEGRAL OF
C                       PLANCK FUNCTION FROM ZERO TO WNUMLO (I=1) OR
C                       WNUMHI (I=2)
C        PI       :  3.14159...
C        SIGMA    :  STEFAN-BOLTZMANN CONSTANT (W/M**2/K**4)
C        SIGDPI   :  SIGMA / PI
C        SMALLV   :  NUMBER OF TIMES THE POWER SERIES IS USED (0,1,2)
C        V(I)     :  C2 * (WNUMLO(I=1) OR WNUMHI(I=2)) / TEMPERATURE
C        VCUT     :  POWER-SERIES CUTOFF POINT
C        VCP      :  EXPONENTIAL SERIES CUTOFF POINTS
C        VMAX     :  LARGEST ALLOWABLE ARGUMENT OF 'EXP' FUNCTION
C
      PARAMETER  ( A1 = 1./3., A2 = -1./8., A3 = 1./60., A4 = -1./5040.,
     $             A5 = 1./272160., A6 = -1./13305600. )
      DOUBLE PRECISION d1mach
      INTEGER  SMALLV
      REAL     C2, CONC, D(2), EPSIL, EX, MV, P(2), SIGMA, SIGDPI,
     $         V(2), VCUT, VCP(7), VSQ
      SAVE     CONC, VMAX, EPSIL, SIGDPI
      DATA     c1 /1.1911E-8/
      DATA     C2 / 1.438786 /,  SIGMA / 5.67032E-8 /,
     $         VCUT / 1.5 /, VCP / 10.25, 5.7, 3.9, 2.9, 2.3, 1.9, 0.0 /
      DATA     PI / 0.0 /
      F(X) = X**3 / ( EXP(X) - 1 )
C
C
      IF ( PI.EQ.0.0 )  THEN
         PI = 2. * ASIN( 1.0 )
         VMAX = ALOG( R1MACH(2) )
         EPSIL = R1MACH(4)
         SIGDPI = SIGMA / PI
         CONC = 15. / PI**4
      END IF
C
      IF( T.LT.0.0 .OR. WNUMHI.LT.WNUMLO .OR. WNUMLO.LT.0. ) 
     $   CALL ERRMSG( 'PLKAVG--TEMPERATURE OR WAVENUMS. WRONG', .TRUE.)
C
      IF ( T.LT.1.E-4 )  THEN
         TPLKAVG = 0.0
         RETURN
      ENDIF
C
*
      IF ( wnumhi .eq. wnumlo ) THEN
         wvn  =  wnumhi
         arg  = EXP( - C2 * wvn / T )
         tplkavg = c1 * (wvn**3.) * arg / ( 1. - arg )
         RETURN
      ENDIF
*     
      V(1) = C2 * WNUMLO / T
      V(2) = C2 * WNUMHI / T
      IF ( V(1).GT.EPSIL .AND. V(2).LT.VMAX .AND.
     $     (WNUMHI-WNUMLO)/WNUMHI .LT. 1.E-2 )  THEN
C
C                          ** WAVENUMBERS ARE VERY CLOSE.  GET INTEGRAL
C                          ** BY ITERATING SIMPSON RULE TO CONVERGENCE.
         HH = V(2) - V(1)
         OLDVAL = 0.0
         VAL0 = F( V(1) ) + F( V(2) )
C
         DO  2  N = 1, 10
            DEL = HH / (2*N)
            VAL = VAL0
            DO  1  K = 1, 2*N-1
               VAL = VAL + 2*(1+MOD(K,2)) * F( V(1) + K*DEL )
    1       CONTINUE
            VAL = DEL/3. * VAL
            IF ( ABS( (VAL-OLDVAL)/VAL ) .LE. 1.E-6 )  GO TO 3
            OLDVAL = VAL
    2    CONTINUE
         CALL ERRMSG( 'PLKAVG--SIMPSON RULE DIDNT CONVERGE', .FALSE. )
C
    3    TPLKAVG = SIGDPI * T**4 * CONC * VAL
         RETURN
      END IF
C
      SMALLV = 0
      DO  50  I = 1, 2
C
         IF( V(I).LT.VCUT )  THEN
C                                   ** USE POWER SERIES
            SMALLV = SMALLV + 1
            VSQ = V(I)**2
            P(I) =  CONC * VSQ * V(I) * ( A1 + V(I) * ( A2 + V(I) *
     $                ( A3 + VSQ * ( A4 + VSQ * ( A5 + VSQ*A6 ) ) ) ) )
         ELSE
C                    ** USE EXPONENTIAL SERIES
            MMAX = 0
C                                ** FIND UPPER LIMIT OF SERIES
   20       MMAX = MMAX + 1
               IF ( V(I).LT.VCP( MMAX ) )  GO TO 20
C
            EX = EXP( - V(I) )
            EXM = 1.0
            D(I) = 0.0
C
            DO  30  M = 1, MMAX
               MV = M * V(I)
               EXM = EX * EXM
               D(I) = D(I) +
     $                EXM * ( 6. + MV*( 6. + MV*( 3. + MV ) ) ) / M**4
   30       CONTINUE
C
            D(I) = CONC * D(I)
         END IF
C
   50 CONTINUE
C
      IF ( SMALLV .EQ. 2 ) THEN
C                                    ** WNUMLO AND WNUMHI BOTH SMALL
         TPLKAVG = P(2) - P(1)
C
      ELSE IF ( SMALLV .EQ. 1 ) THEN
C                                    ** WNUMLO SMALL, WNUMHI LARGE
         TPLKAVG = 1. - P(1) - D(2)
C
      ELSE
C                                    ** WNUMLO AND WNUMHI BOTH LARGE
         TPLKAVG = D(1) - D(2)
C
      END IF
C
      TPLKAVG = SIGDPI * T**4 * TPLKAVG
      IF( TPLKAVG.EQ.0.0 )
     $    CALL ERRMSG( 'PLKAVG--RETURNS ZERO; POSSIBLE UNDERFLOW',
     $                 .FALSE. )
C
      RETURN
      END
      SUBROUTINE tprtinp( albedo, btemp, deltam, dtauc, fbeam, fisot,   
     $     flyr, gg, lyrcut, nlyr, planck, nstr, ntau, oprim,
     $     spher, ssalb, tauc, taucpr, temis, temper, ttemp, umu0,
     $     utau, wvnmlo, wvnmhi )
*
* Print values of input variables
*
      LOGICAL  deltam, lyrcut, planck, spher
      REAL     dtauc(*), flyr(*), gg(*), oprim(*), ssalb(*),
     $         tauc( 0:* ), taucpr( 0:* ), temper( 0:* ), utau(*)
*
      WRITE( *,1010 )  nstr, nlyr
      WRITE( *,1030 )  ntau, (utau(lu), lu = 1, ntau)
      IF ( spher ) WRITE(*,1090)
      IF ( .NOT. planck  )  WRITE( *,1100 )
      WRITE( *,1060 ) fbeam, umu0,  fisot
      WRITE( *,1080 ) albedo
      IF ( planck )  WRITE( *,1110 ) wvnmlo, wvnmhi, btemp,
     $                                    ttemp, temis
      IF ( deltam )       WRITE( *,1120 )
      IF ( .NOT.deltam )  WRITE( *,1130 )
      IF ( lyrcut )       WRITE( *,1170 )
      IF ( planck )       WRITE ( *,1190 )
      IF ( .NOT. planck ) WRITE ( *,1191 )
      yessct = 0.0
      DO 10 lc = 1, nlyr
         yessct = yessct + ssalb(lc)
         IF( planck )
     $       WRITE( *,1200 )  lc, dtauc(lc), tauc(lc), ssalb(lc),
     $        flyr(lc), taucpr(lc)-taucpr(lc-1), taucpr(lc),
     $        oprim(lc), gg(lc), temper(lc-1)
         IF( .NOT.  planck )
     $       WRITE( *,1200 )  lc, dtauc(lc), tauc(lc), ssalb(lc),
     $        flyr(lc), taucpr(lc)-taucpr(lc-1), taucpr(lc),
     $        oprim(lc), gg(lc)
 10   CONTINUE
      IF( planck )  WRITE( *,1210 ) temper(nlyr)
*
      RETURN
*
1010  FORMAT ( /, ' No. streams =', I4,
     $  '     No. computational layers =', I4 )
1030  FORMAT( I4,' User optical depths :',10F10.4, /, (26X,10F10.4) )
1060  FORMAT( '    Incident beam with intensity =', 1P,E11.3, ' and',
     $ ' polar angle cosine = ', 0P,F8.5,
     $ /,'    plus isotropic incident intensity =', 1P,E11.3 )
1080  FORMAT( '    Bottom albedo (lambertian) =', 0P,F8.4 )
1090  FORMAT( ' Pseudo-spherical geometry invoked' )
1100  FORMAT( ' No thermal emission' )
1110  FORMAT( '    Thermal emission in wavenumber interval :', 2F14.4,/,
     $   '    bottom temperature =', F10.2, '     top temperature =',
     $   F10.2,'    top emissivity =', F8.4 )
1120  FORMAT( ' Uses delta-m method' )
1130  FORMAT( ' Does not use delta-m method' )
1150  FORMAT( ' Calculate fluxes and intensities' )
1170  FORMAT( ' Sets radiation = 0 below absorption optical depth 10' )
1190  FORMAT( /, 37X, '<------------- delta-m --------------->', /,
     $'                   total    single                           ',
     $               'total    single', /,
     $'       optical   optical   scatter   truncated   ',
     $   'optical   optical   scatter    asymm', /,
     $'         depth     depth    albedo    fraction     ',
     $     'depth     depth    albedo   factor   temperature' )
1191  FORMAT( /, 37X, '<------------- delta-m --------------->', /,
     $'                   total    single                           ',
     $               'total    single', /,
     $'       optical   optical   scatter   truncated   ',
     $   'optical   optical   scatter    asymm', /,
     $'         depth     depth    albedo    fraction     ',
     $     'depth     depth    albedo   factor' )
1200  FORMAT( I4, 2F10.4, F10.5, F12.5, 2F10.4, F10.5, F9.4,F14.3 )
1210  FORMAT( 85X, F14.3 )
1300  FORMAT( I6, 10F11.6, /, (6X,10F11.6) )
*
      END
*
      SUBROUTINE twostr( albedo, btemp, deltam, dtauc, fbeam, fisot, 
     $     gg, header, ierror, maxcly, maxulv, newgeo, nlyr, planck,
     $     ntau, prnt, quiet, radius, spher, ssalb, temis,  
     $     temper, ttemp, umu0,  usrtau, utau, wvnmlo,
     $     wvnmhi, zd, dfdt, flup, rfldir, rfldn, uavg )
************************************************************************
* Copyright (C) 1993, 1994, 1995 Arve Kylling
*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 1, or (at your option)
* any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY of FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* To obtain a copy of the GNU General Public License write to the
* Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139,
* USA.
*
************************************************************************
*     twostr solves the radiative transfer equation in the two-stream
*     approximation. Based on the general-purpose algorithm
*     DISORT, but both simplified and extended as follows:
*
*     1) The only quantities calculated are mean intensities
*        and fluxes.
*         
*     2) The medium may be taken to be pseudo-spherical (spher=.TRUE.).
*
*     3) Only Lambertian reflection at the bottom boundary is allowed
*
*    ( See twostr.doc for complete documentation )
*
************************************************************************
*
*                 I/O variable specifications
*
      CHARACTER  header*127
      LOGICAL deltam, newgeo, planck, prnt(2), quiet, spher, usrtau
      INTEGER ierror(22), maxcly, maxulv, nlyr, ntau
      REAL    albedo, btemp, dtauc( maxcly ), fbeam, fisot, gg(maxcly),
     $        radius, ssalb( maxcly ), temper( 0:maxcly ), temis, ttemp,
     $        wvnmlo, wvnmhi, umu0, utau( maxulv ), zd( 0:maxcly )
*
      REAL    rfldir( maxulv ), rfldn( maxulv ), flup( maxulv ),
     $        uavg( maxulv ), dfdt( maxulv )
*
*+---------------------------------------------------------------------+
*      Routines called (in order):  tzeroal, tchekin, settwo, tprtinp,
*                                   hopsol, setmtx, solbcd, tfluxes,
*+---------------------------------------------------------------------+
*
*  Index conventions (for all do-loops and all variable descriptions):
*
*     iq     :  For quadrature angles
*
*     lu     :  For user levels
*
*     lc     :  For computational layers (each having a different
*               single-scatter albedo and/or phase function)
*
*     lev    :  For computational levels
*
*     ls     :  Runs from 0 to 2*mxcly+1, ls = 1,2,3 refers to top,
*               center and bottom of layer 1, ls = 3,4,5 refers to
*               top, center and bottom of layer 2, etc.....
*
*+---------------------------------------------------------------------+
*               I n t e r n a l    v a r i a b l e s
*
*   b()               Right-hand side vector of Eqs. KST(38-41), set 
*                     in *solbcd*
*
*   bplank            Intensity emitted from bottom boundary
*
*   cband()           Matrix of left-hand side of the linear system
*                     Eqs. KST(38-41);  in tridiagonal form
*
*   ch(lc)            The Chapman-factor to correct for pseudo-
*                     spherical geometry in the direct beam.
*
*   chtau(lc)         The optical depth in spherical geometry.
*
*   cmu               Computational polar angle, single or double
*                     Gaussian quadrature rule used, see routine
*                     -settwo-
*
*   expbea(lc)        Transmission of direct beam in delta-m optical
*                     depth coordinates
*
*   fldn(lu)          Diffuse down flux (delta-m scaled)
*
*   fldir(lu)         Direct beam flux (delta-m scaled)
*
*   flyr(lc)          Truncated fraction in delta-m method
*
*   kk(lc)            Eigenvalues in Eq. KST(20)
*
*   layru(lu)         Computational layer in which user output level
*                     -utau(lu)- is located
*
*   ll(iq,lc)         Constants of integration C-tilde in Eqs. KST(42-43)
*                     obtained by solving Eqs. KST(38-41)
*
*   lyrcut            True, radiation is assumed zero below layer
*                     -ncut- because of almost complete absorption
*
*   ncut              Computational layer number in which absorption
*                     optical depth first exceeds -abscut-
*
*   oprim(lc)         Single scattering albedo after delta-m scaling
*
*   pass1             True on first entry, false thereafter
*
*   pkag(0:lc)        Integrated Planck function for internal emission
*                     at layer boundaries
*
*   pkagc(lc)         Integrated Planck function for internal emission
*                     at layer center
*
*   rr(lc)            Eigenvectors at polar quadrature angles.
*
*   tauc(0:lc)        Cumulative optical depth (un-delta-m-scaled)
*
*   taucpr(0:lc)      Cumulative optical depth (delta-m-scaled if
*                     deltam = true, otherwise equal to -tauc-)
*
*   tplank            Intensity emitted from top boundary
*
*   u0c(iq,lu)        Azimuthally-averaged intensity
*
*   utaupr(lu)        Optical depths of user output levels in delta-m
*                     coordinates;  equal to  -utau(lu)- if no delta-m
*
*   xb_0d(lc)         x-sub-zero-sup-minus in expansion of pseudo 
*                     spherical beam source, Eq. KST(22)
*
*   xb_0u(lc)         x-sub-zero-sup-plus in expansion of pseudo 
*                     spherical beam source, Eq. KST(22)
*
*   xb_1d(lc)         x-sub-one-sup-minus in expansion of pseudo 
*                     spherical beam source, Eq. KST(22)
*
*   xb_1u(lc)         x-sub-one-sup-plus in expansion of pseudo 
*                     spherical beam source, Eq. KST(22)
*
*   xp_0(lc)          x-sub-zero in expansion of thermal source func-
*                     tion; see Eq. KST(22) (has no (mu) dependence)
*
*   xp_1(lc)          x-sub-one in expansion of thermal source func-
*                     tion; see Eq. KST(22) (has no (mu) dependence)
*
*   yb_0d(lc)         y-sub-zero-sup-minus in Eq. KST(23), solution for
*                     pseudo spherical beam source.
*
*   yb_0u(lc)         y-sub-zero-sup-plus in Eq. KST(23), solution for
*                     pseudo spherical beam source.
*
*   yb_1d(lc)         y-sub-one-sup-minus in Eq. KST(23), solution for
*                     pseudo spherical beam source.
*
*   yb_1u(lc)         y-sub-one-sup-plus in Eq. KST(23), solution for
*                     pseudo spherical beam source.
*
*   yp_0d(lc)         y-sub-zero-sup-minus in Eq. KST(23), solution for
*                     thermal source.
*
*   yp_0u(lc)         y-sub-zero-sup-plus in Eq. KST(23), solution for
*                     thermal source.
*
*   yp_1d(lc)         y-sub-one-sup-minus in Eq. KST(23), solution for
*                     thermal source.
*
*   yp_1u(lc)         y-sub-one-sup-plus in Eq. KST(23), solution for
*                     thermal source.
*
*   zb_a(lc)          Alfa coefficient in Eq.  KST(22) for pseudo-
*                     spherical beam source.
*
*   zp_a(lc)          Alfa coefficient in Eq. KST(22) for thermal
*                     source.
*
c +-------------------------------------------------------------------+
c
c  LOCAL SYMBOLIC DIMENSIONS 
c
c  MXCLY  = Max no. of computational layers
c  MXULV  = Max no. of output levels
c  MXCMU  = Max no. of computation polar angles
c  MXUMU  = Max no. of output polar angles
c  MXPHI  = Max no. of output azimuthal angles
c +-------------------------------------------------------------------+

      INTEGER MXCLY, MXULV, MXTCMU, MI9M2, NNLYRI
      INCLUDE "DISORT.MXD"
      PARAMETER ( MXTCMU = 2, MI9M2 = 4,
     &          NNLYRI = MXTCMU*MXCLY)
*
      LOGICAL lyrcut, pass1
      INTEGER ipvt(nnlyri), iret, layru(mxulv), nerr
      REAL b(nnlyri), cband(mi9m2,nnlyri), ch(mxcly),
     $     chtau(0:2*mxcly+1), cmu, dtaucpr(mxcly), expbea(0:mxcly),
     $     flyr(mxcly), fldn(mxulv), fldir(mxulv), ggprim(mxcly),
     $     kk(mxcly),  ll(mxtcmu,mxcly), oprim( mxcly ), pkag(0:mxcly),
     $     pkagc(mxcly), rr(mxcly), tauc(0:mxcly),
     $     taucpr(0:mxcly), u0c(mxtcmu,mxulv), utaupr(mxulv), 
     $     xb_0d(mxcly), xb_0u(mxcly), xb_1d(mxcly), xb_1u(mxcly), 
     $     xp_0(mxcly), xp_1(mxcly), yb_0d(mxcly), yb_0u(mxcly), 
     $     yb_1d(mxcly), yb_1u(mxcly), yp_0d(mxcly), yp_0u(mxcly),
     $     yp_1d(mxcly), yp_1u(mxcly), zb_a(mxcly),
     $     zp_a(mxcly), subd(2*mxcly),
     $     diag(2*mxcly), superd(2*mxcly), sol(2*mxcly)
*
      SAVE  cmu, pass1, pi, epsil, biggest, smallest
      DATA  pass1 / .TRUE. /, nerr / 22 /
*
      IF ( pass1 )  THEN
         pi = 2. * ASIN(1.0)
         epsil =  0.00001             ! Typical value for 32 bits machine
         pass1 = .FALSE.
      END IF
*
      IF ( prnt(1) )  WRITE( *,1010 )  header
*
* Zero some arrays (not strictly necessary, but otherwise unused 
* parts of arrays collect garbage)
*
      DO i = 1, nerr
         ierror(i) = 0
      ENDDO
      CALL  zeroit( ch    , mxcly )
      CALL tzeroal( kk, ll, mxtcmu, mxcly,
     $     rr, xb_0d, xb_0u, xb_1d, xb_1u,
     $     xp_0, xp_1, yb_0d, yb_0u, yb_1d, yb_1u, yp_0d,
     $     yp_0u, yp_1d, yp_1u, zb_a, zp_a )
*
* Calculate cumulative optical depth and dither single-scatter
* albedo to improve numerical behavior of eigenvalue/vector
* computation
*
      tauc( 0 ) = 0.
      CALL  zeroit( tauc(0), mxcly+1 )
      DO 20  lc = 1, nlyr
         IF( ssalb(lc) .GT. 1.0 - epsil )  ssalb(lc) = 1.0 - epsil
         tauc(lc) = tauc(lc-1) + dtauc(lc)
 20   CONTINUE
*
* Check input dimensions and variables
*
      CALL tchekin( albedo, btemp, dtauc, fbeam, fisot, gg,
     $     ierror, maxcly, maxulv, mxcly, mxulv, nlyr, ntau, planck,
     $     quiet, spher, ssalb, tauc, temis, temper, ttemp, 
     $     umu0, usrtau, utau, wvnmlo, wvnmhi, zd )
*
      iret = 0
      DO ierr = 1, nerr
         IF ( ierror(ierr) .NE. 0 ) THEN
            iret = 1
            WRITE(*,'(/,A,I4,/)')  "TWOSTR REPORTS FATAL ERROR: ",
     $           ierr
         ENDIF
      ENDDO
      IF ( iret .EQ. 1 ) RETURN
*     
* Perform various setup operations
*
      CALL settwo(biggest, bplank, btemp, ch, chtau, cmu, deltam,  
     $     dtauc, dtaucpr, expbea, fbeam, flyr, gg, ggprim,
     $     layru, lyrcut, ncut, newgeo, nlyr, nn, nstr, ntau, 
     $     oprim, pi, pkag, pkagc, planck, radius, smallest, spher, 
     $     ssalb, tauc, taucpr, temis, temper, tplank, ttemp, umu0, 
     $     usrtau, utau, utaupr, wvnmlo, wvnmhi, zd )
*     
* Print input information
*
      IF ( prnt(1) )
     $     CALL tprtinp( albedo, btemp, deltam, dtauc, fbeam, fisot,   
     $     flyr, gg, lyrcut, nlyr, planck, nstr, ntau, oprim,
     $     spher, ssalb, tauc, taucpr, temis, temper, ttemp, umu0,
     $     utau, wvnmlo, wvnmhi )
*
* Calculate the homogenous and particular solutions
*
      CALL hopsol( biggest, ch, chtau, cmu, fbeam, ggprim, kk, ncut,
     $     oprim, pi, pkag, pkagc, planck, rr, smallest,
     $     taucpr, umu0, xb_0d, xb_0u, xb_1d, xb_1u, xp_0,xp_1,
     $     yb_0d, yb_0u, yb_1d, yb_1u, yp_0d, yp_0u, yp_1d,yp_1u, zb_a,
     $     zp_a )
*     
* Solve for constants of integration in homogeneous solution
* (general boundary conditions)
*
      CALL solbcd( albedo, b, bplank, cband, cmu, diag, expbea,
     $     fbeam, fisot, ipvt, kk, ll, lyrcut, mi9m2, mxtcmu,
     $     ncut, nnlyri, pi, rr, subd,
     $     superd, sol, tplank, taucpr, umu0, yb_0d, yb_0u, yb_1d,
     $     yb_1u, yp_0d, yp_0u, yp_1d, yp_1u, zb_a,zp_a )
*
* Compute upward and downward fluxes, mean intensities and
* flux divergences.
*
      CALL tfluxes( ch, cmu, dfdt, fbeam, flup, fldn, fldir, kk,
     $     layru, ll, lyrcut, maxulv, mxtcmu, mxulv, ncut, ntau,
     $     pi, planck, prnt, oprim, rfldir, rfldn, rr, spher, ssalb,
     $     taucpr, u0c, uavg, umu0, utau, utaupr, xp_0, xp_1,
     $     yb_0d, yb_0u, yb_1d, yb_1u, yp_0d, yp_0u, yp_1d, yp_1u,
     $     zb_a, zp_a )
*     
      RETURN
*
1010  FORMAT ( ////, 1X, 120('*'), /, 25X,
     $  'Two stream method radiative transfer program, version 1.20',
     $  /, 1X, A, /, 1X, 120('*') )
*
      END
*
      SUBROUTINE tzeroal( kk, ll, mxtcmu, mxcly,
     $     rr, xb_0d, xb_0u, xb_1d, xb_1u,
     $     xp_0, xp_1, yb_0d, yb_0u, yb_1d, yb_1u, yp_0d,
     $     yp_0u, yp_1d, yp_1u, zb_a, zp_a )
*
* Zero arrays
*
      REAL  kk(*), ll(mxtcmu,*), rr(*), 
     $      xb_0d(*), xb_0u(*), xb_1d(*),
     $      xb_1u(*), xp_0(*), xp_1(*), yb_0d(*), yb_0u(*),
     $      yb_1d(*), yb_1u(*), yp_0d(*), yp_0u(*), yp_1d(*), 
     $      yp_1u(*), zb_a(*), zp_a(*)
*
      CALL  zeroit( kk,     mxcly )
      CALL  zeroit( ll,     mxtcmu*mxcly )
      CALL  zeroit( rr,     mxcly )
      CALL  zeroit( xb_0d, mxcly )
      CALL  zeroit( xb_0u, mxcly )
      CALL  zeroit( xb_1d, mxcly )
      CALL  zeroit( xb_1u, mxcly )
      CALL  zeroit( xp_0, mxcly )
      CALL  zeroit( xp_1, mxcly )
      CALL  zeroit( yb_0d, mxcly )
      CALL  zeroit( yb_0u, mxcly )
      CALL  zeroit( yb_1d, mxcly )
      CALL  zeroit( yb_1u, mxcly )
      CALL  zeroit( yp_0d, mxcly )
      CALL  zeroit( yp_0u, mxcly )
      CALL  zeroit( yp_1d, mxcly )
      CALL  zeroit( yp_1u, mxcly )
      CALL  zeroit( zb_a, mxcly )
      CALL  zeroit( zp_a, mxcly )
*
      RETURN
      END
*

