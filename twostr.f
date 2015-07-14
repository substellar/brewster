      SUBROUTINE twostr( albedo, btemp, deltam, dtauc, fbeam, fisot, 
     $     gg, header, ierror, maxcly, maxulv, nlyr, planck,
     $     ntau, prnt, quiet, radius, spher, ssalb, temis,  
     $     temper, ttemp, umu0,  usrtau, utau, wvnmlo,
     $     wvnmhi, zd, dfdt, flup, rfldir, rfldn, uavg )
c***********************************************************************
c Copyright (C) 1993, 1994, 1995 Arve Kylling
c
c This program is free software; you can redistribute it and/or modify
c it under the terms of the GNU General Public License as published by
c the Free Software Foundation; either version 1, or (at your option)
c any later version.
c
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY of FITNESS FOR A PARTICULAR PURPOSE. See the
c GNU General Public License for more details.
c
c To obtain a copy of the GNU General Public License write to the
c Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139,
c USA.
c
c***********************************************************************
c     twostr solves the radiative transfer equation in the two-stream
c     approximation. Based on the general-purpose algorithm
c     DISORT, but both simplified and extended as follows:
c
c     1) The only quantities calculated are mean intensities
c        and fluxes.
c         
c     2) The medium may be taken to be pseudo-spherical (spher=.TRUE.).
c
c     3) Only Lambertian reflection at the bottom boundary is allowed
c
c    ( See twostr.doc for complete documentation )
c
c***********************************************************************
c
c                 I/O variable specifications
c
      CHARACTER  header*127
      LOGICAL deltam, planck, prnt(2), quiet, spher, usrtau
      INTEGER ierror(22), maxcly, maxulv, nlyr, ntau
      REAL    albedo, btemp, dtauc( maxcly ), fbeam, fisot, gg(maxcly),
     $        radius, ssalb( maxcly ), temper( 0:maxcly ), temis, ttemp,
     $        wvnmlo, wvnmhi, umu0, utau( maxulv ), zd( 0:maxcly )
c
      REAL    rfldir( maxulv ), rfldn( maxulv ), flup( maxulv ),
     $        uavg( maxulv ), dfdt( maxulv )
c
c+---------------------------------------------------------------------+
c      Routines called (in order):  tzeroal, tchekin, settwo, tprtinp,
c                                   hopsol, setmtx, solbcd, tfluxes,
c+---------------------------------------------------------------------+
c
c  Index conventions (for all do-loops and all variable descriptions):
c
c     iq     :  For quadrature angles
c
c     lu     :  For user levels
c
c     lc     :  For computational layers (each having a different
c               single-scatter albedo and/or phase function)
c
c     lev    :  For computational levels
c
c     ls     :  Runs from 0 to 2*mxcly+1, ls = 1,2,3 refers to top,
c               center and bottom of layer 1, ls = 3,4,5 refers to
c               top, center and bottom of layer 2, etc.....
c
c+---------------------------------------------------------------------+
c               I n t e r n a l    v a r i a b l e s
c
c   b()               Right-hand side vector of Eqs. KST(38-41), set 
c                     in *solbcd*
c
c   bplank            Intensity emitted from bottom boundary
c
c   cband()           Matrix of left-hand side of the linear system
c                     Eqs. KST(38-41);  in tridiagonal form
c
c   ch(lc)            The Chapman-factor to correct for pseudo-
c                     spherical geometry in the direct beam.
c
c   chtau(lc)         The optical depth in spherical geometry.
c
c   cmu               Computational polar angle, single or double
c                     Gaussian quadrature rule used, see routine
c                     -settwo-
c
c   expbea(lc)        Transmission of direct beam in delta-m optical
c                     depth coordinates
c
c   fldn(lu)          Diffuse down flux (delta-m scaled)
c
c   fldir(lu)         Direct beam flux (delta-m scaled)
c
c   flyr(lc)          Truncated fraction in delta-m method
c
c   kk(lc)            Eigenvalues in Eq. KST(20)
c
c   layru(lu)         Computational layer in which user output level
c                     -utau(lu)- is located
c
c   ll(iq,lc)         Constants of integration C-tilde in Eqs. KST(42-43)
c                     obtained by solving Eqs. KST(38-41)
c
c   lyrcut            True, radiation is assumed zero below layer
c                     -ncut- because of almost complete absorption
c
c   ncut              Computational layer number in which absorption
c                     optical depth first exceeds -abscut-
c
c   oprim(lc)         Single scattering albedo after delta-m scaling
c
c   pass1             True on first entry, false thereafter
c
c   pkag(0:lc)        Integrated Planck function for internal emission
c                     at layer boundaries
c
c   pkagc(lc)         Integrated Planck function for internal emission
c                     at layer center
c
c   rr(lc)            Eigenvectors at polar quadrature angles.
c
c   tauc(0:lc)        Cumulative optical depth (un-delta-m-scaled)
c
c   taucpr(0:lc)      Cumulative optical depth (delta-m-scaled if
c                     deltam = true, otherwise equal to -tauc-)
c
c   tplank            Intensity emitted from top boundary
c
c   u0c(iq,lu)        Azimuthally-averaged intensity
c
c   utaupr(lu)        Optical depths of user output levels in delta-m
c                     coordinates;  equal to  -utau(lu)- if no delta-m
c
c   xb_0d(lc)         x-sub-zero-sup-minus in expansion of pseudo 
c                     spherical beam source, Eq. KST(22)
c
c   xb_0u(lc)         x-sub-zero-sup-plus in expansion of pseudo 
c                     spherical beam source, Eq. KST(22)
c
c   xb_1d(lc)         x-sub-one-sup-minus in expansion of pseudo 
c                     spherical beam source, Eq. KST(22)
c
c   xb_1u(lc)         x-sub-one-sup-plus in expansion of pseudo 
c                     spherical beam source, Eq. KST(22)
c
c   xp_0(lc)          x-sub-zero in expansion of thermal source func-
c                     tion; see Eq. KST(22) (has no (mu) dependence)
c
c   xp_1(lc)          x-sub-one in expansion of thermal source func-
c                     tion; see Eq. KST(22) (has no (mu) dependence)
c
c   yb_0d(lc)         y-sub-zero-sup-minus in Eq. KST(23), solution for
c                     pseudo spherical beam source.
c
c   yb_0u(lc)         y-sub-zero-sup-plus in Eq. KST(23), solution for
c                     pseudo spherical beam source.
c
c   yb_1d(lc)         y-sub-one-sup-minus in Eq. KST(23), solution for
c                     pseudo spherical beam source.
c
c   yb_1u(lc)         y-sub-one-sup-plus in Eq. KST(23), solution for
c                     pseudo spherical beam source.
c
c   yp_0d(lc)         y-sub-zero-sup-minus in Eq. KST(23), solution for
c                     thermal source.
c
c   yp_0u(lc)         y-sub-zero-sup-plus in Eq. KST(23), solution for
c                     thermal source.
c
c   yp_1d(lc)         y-sub-one-sup-minus in Eq. KST(23), solution for
c                     thermal source.
c
c   yp_1u(lc)         y-sub-one-sup-plus in Eq. KST(23), solution for
c                     thermal source.
c
c   zb_a(lc)          Alfa coefficient in Eq.  KST(22) for pseudo-
c                     spherical beam source.
c
c   zp_a(lc)          Alfa coefficient in Eq. KST(22) for thermal
c                     source.
c
c+---------------------------------------------------------------------+
c   Local symbolic dimensions:
c
c       mxcly  = max no. of computational layers
c       mxulv  = max no. of output levels
c       mxcmu  = max no. of computation polar angles (=2 in two-stream)
c+---------------------------------------------------------------------+
      PARAMETER ( mxcly = 61, mxulv = 1, mxcmu = 2, mi = mxcmu/2, 
     $            mi9m2 = 4, nnlyri = mxcly*mxcmu )
c
      LOGICAL lyrcut, pass1
      INTEGER ipvt(nnlyri), iret, layru(mxulv), nerr
      REAL b(nnlyri), cband(mi9m2,nnlyri), ch(mxcly), chtau(2*mxcly+1),
     $     cmu, dtaucpr(mxcly), expbea(0:mxcly),
     $     flyr(mxcly), fldn(mxulv), fldir(mxulv), ggprim(mxcly),
     $     kk(mxcly),  ll(mxcmu,mxcly), oprim( mxcly ), pkag(0:mxcly),
     $     pkagc(mxcly), rr(mxcly), tauc(0:mxcly),
     $     taucpr(0:mxcly), u0c(mxcmu,mxulv), utaupr(mxulv), 
     $     xb_0d(mxcly), xb_0u(mxcly), xb_1d(mxcly), xb_1u(mxcly), 
     $     xp_0(mxcly), xp_1(mxcly), yb_0d(mxcly), yb_0u(mxcly), 
     $     yb_1d(mxcly), yb_1u(mxcly), yp_0d(mxcly), yp_0u(mxcly),
     $     yp_1d(mxcly), yp_1u(mxcly), zb_a(mxcly),
     $     zp_a(mxcly), subd(2*mxcly),
     $     diag(2*mxcly), superd(2*mxcly), sol(2*mxcly) 
c
      SAVE  pass1, pi, epsil
      DATA  pass1 / .TRUE. /, nerr / 22 /
c
      IF ( pass1 )  THEN
         pi = 2. * ASIN(1.0)
         epsil =  0.00001             ! Typical value for 32 bits machine
         pass1 = .FALSE.
      END IF
c
      IF ( prnt(1) )  WRITE( *,1010 )  header
c
c Zero some arrays (not strictly necessary, but otherwise unused 
c parts of arrays collect garbage)
c
      DO i = 1, nerr
         ierror(i) = 0
      ENDDO
      CALL tzeroal( ch, kk, ll, mxcmu, mxcly,
     $     nnlyri, rr, xb_0d, xb_0u, xb_1d, xb_1u,
     $     xp_0, xp_1, yb_0d, yb_0u, yb_1d, yb_1u, yp_0d,
     $     yp_0u, yp_1d, yp_1u, zb_a, zp_a )
c
c Calculate cumulative optical depth and dither single-scatter
c albedo to improve numerical behavior of eigenvalue/vector
c computation
c
      tauc( 0 ) = 0.
      CALL  tzeroit( tauc(0), mxcly+1 )
      DO 20  lc = 1, nlyr
         IF( ssalb(lc) .GT. 1.0 - epsil )  ssalb(lc) = 1.0 - epsil
         tauc(lc) = tauc(lc-1) + dtauc(lc)
 20   CONTINUE
c
c Check input dimensions and variables
c
      CALL tchekin( albedo, btemp, dtauc, fbeam, fisot, gg,
     $     ierror, maxcly, maxulv, mxcly, mxulv, nlyr, ntau, planck,
     $     quiet, spher, ssalb, tauc, temis, temper, ttemp, 
     $     umu0, usrtau, utau, wvnmlo, wvnmhi, zd )
c
      iret = 0
      DO ierr = 1, nerr
         IF ( ierror(ierr) .NE. 0 ) THEN
            iret = 1
            IF ( .NOT. quiet )  THEN
               WRITE(*,'(/,A,I4,/)')  "TWOSTR REPORTS FATAL ERROR: ",
     $              ierr
            ENDIF
         ENDIF
      ENDDO
      IF ( iret .EQ. 1 ) RETURN
c     
c Perform various setup operations
c
      CALL settwo( albedo, biggest, bplank, btemp, ch, chtau, cmu,
     $     deltam,  
     $     dtauc, dtaucpr, expbea, fbeam, flyr, gg, ggprim,
     $     layru, lyrcut, mxcly, ncut, nlyr, nn, nstr, ntau, oprim,
     $     pkag, pkagc, planck, radius, smallest, spher, ssalb, tauc, 
     $     taucpr, temis, temper, tplank, ttemp, umu0, usrtau,
     $     utau, utaupr, wvnmlo, wvnmhi, zd )
c     
c Print input information
c
      IF ( prnt(1) )
     $     CALL tprtinp( albedo, btemp, deltam, dtauc, fbeam, fisot,   
     $     flyr, gg, lyrcut, nlyr, planck, nstr, ntau, oprim,
     $     spher, ssalb, tauc, taucpr, temis, temper, ttemp, umu0,
     $     utau, wvnmlo, wvnmhi )
c
c Calculate the homogenous and particular solutions
c
      CALL hopsol( biggest, ch, chtau, cmu, fbeam, ggprim, kk, ncut,
     $     nlyr,oprim, pi, pkag, pkagc, planck, radius, rr, smallest,
     $     spher,taucpr, umu0, xb_0d, xb_0u, xb_1d, xb_1u, xp_0,xp_1,
     $     yb_0d, yb_0u, yb_1d, yb_1u, yp_0d, yp_0u, yp_1d,yp_1u, zb_a,
     $     zp_a )
c     
c Solve for constants of integration in homogeneous solution
c (general boundary conditions)
c
      CALL solbcd( albedo, b, bplank, cband, cmu, diag, expbea,
     $     fbeam, fisot, ipvt, kk, ll, lyrcut, mi, mi9m2, mxcmu,
     $     ncol, ncut, nlyr, nn, nnlyri, nstr, pi, rr, subd,
     $     superd, sol, tplank, taucpr, umu0, yb_0d, yb_0u, yb_1d,
     $     yb_1u, yp_0d, yp_0u, yp_1d, yp_1u, zb_a,zp_a )
c
c Compute upward and downward fluxes, mean intensities and
c flux divergences.
c
      CALL tfluxes( ch, cmu, dfdt, fbeam, flup, fldn, fldir, kk,
     $     layru, ll, lyrcut, maxulv, mxcmu, mxulv, ncut, ntau,
     $     pi, planck, prnt, oprim, rfldir, rfldn, rr, spher, ssalb,
     $     taucpr, u0c, uavg, umu0, utau, utaupr, xp_0, xp_1,
     $     yb_0d, yb_0u, yb_1d, yb_1u, yp_0d, yp_0u, yp_1d, yp_1u,
     $     zb_a, zp_a )
c     
      RETURN
c
1010  FORMAT ( ////, 1X, 120('*'), /, 25X,
     $  'Two stream method radiative transfer program, version 1.13',
     $  /, 1X, A, /, 1X, 120('*') )
c
      END
c
      REAL FUNCTION chapmn( lc, taup, tauc, nlyr, zd, dtauc,
     $     zenang, r )
c
c Calculates the Chapman-factor
c
c I n p u t       v a r i a b l e s:
c
c      lc        : Computational layer
c      nlyr      : Number of layers in atmospheric model
c      zd(lc)    : lc = 0, nlyr. zd(lc) is distance from bottom 
c                  surface to top of layer lc. zd(nlyr) = 0.0 km
c      dtauc     : Optical thickness of layer lc (un-delta-m-scaled)
c      zenang    : Solar zenith angle as seen from bottom surface
c      r         : Radius of earth. NOTE: Use the same dimension as zd,
c                  for instance both in km.
c
c O u t p u t      v a r i a b l e s:
c
c      ch        : Chapman-factor. in a pseudo-spherical atmosphere 
c                  replace EXP( -tau/umu0 ) by EXP( -ch(lc) ) in the
c                  beam source in 
c
c I n t e r n a l     v a r i a b l e s:
c
c      dhj       : delta-h-sub-j in Eq. B2 (DS)
c      dsj       : delta-s-sub-j in Eq. B2 (DS)
c      fact      : =1 for first sum in Eq. B2 (DS)
c                  =2 for second sum in Eq. B2 (DS)
c      rj        : r-sub-j in Eq. B1 (DS)
c      rjp1      : r-sub-j+1 in Eq. B1 (DS)
c      xpsinz    : The length of the line OG in Fig. 1, (DS)
c
c
      REAL dtauc(*) ,tauc(0:*), zd(0:*)
c
      pi     = 2.0 * ASIN( 1.0 )
      zenrad = zenang * pi / 180.0
      xp     = r +  zd(lc) + (zd(lc-1) - zd(lc) ) *
     $     ( tauc(lc) - taup ) / dtauc(lc)
      xpsinz = xp * SIN( zenrad )
c
      IF( (zenang.GT.90.0) .AND. (xpsinz.LT.r) ) THEN
        chapmn = 1.0E+20
        RETURN
      END IF
c
c Find index of layer in which the screening height lies
c
      id = lc
      IF( zenang.GT.90.0 ) THEN
         DO 100 j = lc, nlyr
            IF( (xpsinz.LT.( zd(j-1) + r ) ) .AND.
     $                         (xpsinz.GE.( zd(j) + r )) ) id = j
 100     CONTINUE
      END IF
c
      sum = 0.0
      DO 200 j = 1, id
        fact = 1.0
        fact2= 1.0
c
c Include factor of 2 for zenang .gt. 90, second sum in Eq. B2 (DS)
c
        IF( j.GT.lc ) fact = 2.0
        IF(j.EQ.id .AND. id.EQ.lc .AND. zenang.GT.90.0) fact2 = -1.0
c
        rj = r + zd(j-1)
        rjp1 = r + zd(j)
        IF(j.EQ.lc .AND. id.EQ.lc) rjp1 = xp
c
        dhj = zd(j-1) -zd(j)
        IF(id.GT.lc .AND. j.EQ.id) THEN
           dsj = SQRT(rj*rj - xpsinz*xpsinz )
        ELSE
           dsj = SQRT( rj*rj - xpsinz*xpsinz ) -
     $           fact2 * SQRT( rjp1*rjp1 - xpsinz*xpsinz )
        END IF
c
        sum = sum + dtauc(j)*fact* dsj / dhj
c
 200  CONTINUE
c
c Add third term in Eq. B2 (DS)
c
      IF( id.GT.lc ) THEN
        dhj = zd(lc-1) -zd(lc)
        dsj = SQRT( xp*xp - xpsinz*xpsinz ) -
     $        SQRT( (zd(lc)+r)*(zd(lc)+r) - xpsinz*xpsinz )
        sum = sum + dtauc(lc) * dsj / dhj
      END IF
c
      chapmn = sum
      RETURN
      END
c
      SUBROUTINE  tchekin(albedo, btemp, dtauc, fbeam, fisot, gg,
     $     ierror, maxcly, maxulv, mxcly, mxulv, nlyr, ntau, planck,
     $     quiet, spher, ssalb, tauc, temis, temper, ttemp, umu0,
     $     usrtau, utau, wvnmlo, wvnmhi, zd )
c
c Checks the input dimensions and variables
c
      LOGICAL inperr, planck, quiet, spher, usrtau, wrtbad, wrtdim
      INTEGER  maxcly, maxulv, mxcly, mxulv, nlyr, ntau, ierror(*)
      REAL albedo, btemp, dtauc(*), fbeam, fisot, gg(*), 
     $     ssalb(*), tauc(0:*), temis, temper(0:*), ttemp,
     $     umu0, utau(*), wvnmlo, wvnmhi, zd(0:*)
c
      inperr = .FALSE.
      IF ( nlyr.LT.1 ) THEN
         inperr = wrtbad(quiet, 'nlyr' )
         ierror(1) = 1
      ENDIF
      IF ( nlyr.GT.maxcly ) THEN
         inperr = wrtbad(quiet, 'maxcly' )
         ierror(2) = 1
      ENDIF
c
      DO 10  lc = 1, nlyr
         IF ( dtauc(lc).LT.0.0 )  THEN
            inperr = wrtbad( quiet, 'dtauc' )
            ierror(3) = ierror(3) + 1
         ENDIF
         IF ( ssalb(lc).LT.0.0 .OR. ssalb(lc).GT.1.0 ) THEN
            inperr = wrtbad( quiet, 'ssalb' )
            ierror(4) = ierror(4) + 1
         ENDIF
         IF ( planck )  THEN
            IF( lc.EQ.1 .AND. temper(0).LT.0.0 ) THEN
               inperr = wrtbad( quiet, 'temper' )
               ierror(5) = ierror(5) + 1
            ENDIF
            IF( temper(lc).LT.0.0 ) THEN 
               inperr = wrtbad( quiet, 'temper' )
               ierror(5) = ierror(5) + 1
            ENDIF
         ENDIF
         IF( gg(lc).LT.-1.0 .OR. gg(lc).GT.1.0 ) THEN
            inperr = wrtbad( quiet, 'gg' )
            ierror(6) = ierror(6) + 1
         ENDIF
 10   CONTINUE
      IF ( spher ) THEN
         DO 11 lc = 1, nlyr
            IF ( zd(lc) .GT. zd(lc-1) ) THEN
               inperr = wrtbad( quiet, 'zd' )
               ierror(7) = ierror(7) + 1
            ENDIF
 11      CONTINUE  	   
      ENDIF
c
      IF ( usrtau )  THEN
         IF ( ntau.LT.1 ) THEN
            inperr = wrtbad( quiet, 'ntau' )
            ierror(8) = 1
         ENDIF
         IF ( maxulv.LT.ntau ) THEN
            inperr = wrtbad( quiet, 'maxulv' )
            ierror(9) = 1
         ENDIF
         DO 20  lu = 1, ntau
            IF( ABS(utau(lu)-tauc(nlyr)).LE.1.E-4) utau(lu) = tauc(nlyr)
            IF( utau(lu).LT.0.0 .OR. utau(lu).GT.tauc(nlyr) ) THEN
               inperr = wrtbad( quiet, 'utau' )
               ierror(10) = ierror(10) + 1
            ENDIF
 20      CONTINUE
      ELSE
         IF ( maxulv.LT.nlyr+1 ) THEN
            inperr = wrtbad( quiet, 'maxulv' )
            ierror(11) = 1
         ENDIF
      END IF
c
      IF ( fbeam.LT.0.0 ) THEN
         inperr = wrtbad( quiet, 'fbeam' )
         ierror(12) = 1
      ENDIF
      umumin = 0.0
      IF ( spher ) umumin = - 1.0
      IF ( fbeam.GT.0.0 .AND. ( umu0.LE.umumin .OR. umu0.GT.1.0 ) )
     $     THEN
         inperr = wrtbad( quiet, 'umu0' )
         ierror(13) = 1
      ENDIF
      IF ( fisot.LT.0.0 ) THEN
         inperr = wrtbad( quiet, 'fisot' )
         ierror(14) = 1
      ENDIF
      IF ( albedo.LT.0.0 .OR. albedo.GT.1.0 ) THEN
         inperr = wrtbad( quiet, 'albedo' )
         ierror(15) = 1
      ENDIF
c
      IF ( planck )  THEN
         IF ( wvnmlo.LT.0.0 .OR. wvnmhi.LT.wvnmlo ) THEN
            inperr = wrtbad( quiet, 'wvnmlo,hi' )
            ierror(16) = 1
         ENDIF
         IF ( temis.LT.0.0 .OR. temis.GT.1.0 ) THEN
            inperr = wrtbad( quiet, 'temis' )
            ierror(17) = 1
         ENDIF
         IF ( btemp.LT.0.0 ) THEN
            inperr = wrtbad( quiet, 'btemp' )
            ierror(18) = 1
         ENDIF
         IF ( ttemp.LT.0.0 ) THEN
            inperr = wrtbad( quiet, 'ttemp' )
            ierror(19) = 1
         ENDIF
      END IF
c
      IF ( mxcly.LT.nlyr ) THEN
         inperr = wrtdim( quiet, 'mxcly', nlyr )
         ierror(20) = 1
      ENDIF
      IF ( usrtau .AND. mxulv.LT.ntau ) THEN
         inperr = wrtdim( quiet, 'mxulv', ntau )
         ierror(21) = 1
      ENDIF
      IF ( .NOT.usrtau .AND. mxulv.LT.nlyr+1 ) THEN
         inperr = wrtdim( quiet, 'mxulv', nlyr+1 )
         ierror(22) = 1
      ENDIF
c
      IF ( inperr .AND. .NOT. quiet )
     $   CALL errmsg( 'twostr--input and/or dimension errors', .TRUE. )
c
      DO 100  lc = 1, nlyr
         IF ( (planck .AND. ABS(temper(lc)-temper(lc-1)) .GT. 50.0) 
     $          .AND. .NOT. quiet )
     $          CALL errmsg( 'chekin--vertical temperature step may'
     $          // ' be too large for good accuracy', .FALSE. )
100   CONTINUE
c
      RETURN
      END
c
      SUBROUTINE tfluxes( ch, cmu, dfdt, fbeam, flup, fldn, fldir, kk,
     $     layru, ll, lyrcut, maxulv, mxcmu, mxulv, ncut, ntau,
     $     pi, planck, prnt, oprim, rfldir, rfldn, rr, spher, ssalb,
     $     taucpr, u0c, uavg, umu0, utau, utaupr, xp_0, xp_1,
     $     yb_0d, yb_0u, yb_1d, yb_1u, yp_0d, yp_0u, yp_1d, yp_1u,
     $     zb_a, zp_a )
c
c Calculates the radiative fluxes, mean intensity, and flux
c derivative with respect to optical depth from the 
c azimuthally-averaged intensity
c
c I n p u t     v a r i a b l e s:
c
c       ch       :  Chapman factor
c       cmu      :  Abscissae for gauss quadrature over angle cosine
c       kk       :  Eigenvalues 
c       layru    :  Layer number of user level -utau-
c       ll       :  Constants of integration in Eqs. KST(42-43), obtained
c                   by solving Eqs. KST(38-41)
c       lyrcut   :  Logical flag for truncation of comput. layer
c       ncut     :  Number of computational layer where absorption
c                     optical depth exceeds -abscut-
c       oprim    :  Delta-m scaled single scattering albedo
c       rr       :  Eigenvectors at polar quadrature angles 
c       taucpr   :  Cumulative optical depth (delta-m-scaled)
c       utaupr   :  Optical depths of user output levels in delta-m
c                     coordinates;  equal to  -utau- if no delta-m
c       xp_0,    :  Thermal source coefficients x-sub-zero and
c        xp_1         x-sub-one in Eq. KST(22)
c       yb_0d,u, :  Thermal source vectors, Eq. KST(23)
c        yb_1d,u     
c       yp_0d,u, :  Beam source vectors, Eq. KST(23)
c        yp_1d,u     
c       zb_a     :  Beam source coefficient alfa in Eq. KST(22)
c       zp_a     :  Thermal source coefficient alfa in Eq. KST(22)
c       (remainder are 'twostr' input variables)
c
c O u t p u t     v a r i a b l e s:
c
c       u0c      :  Azimuthally averaged intensities at polar
c                     quadrature angle cmu
c       (rfldir, rfldn, flup, dfdt, uavg are 'twostr' output variables)
c
c I n t e r n a l       v a r i a b l e s:
c
c       dirint   :  direct intensity attenuated
c       fdntot   :  total downward flux (direct + diffuse)
c       fldir    :  direct-beam flux (delta-m scaled)
c       fldn     :  diffuse down-flux (delta-m scaled)
c       fnet     :  net flux (total-down - diffuse-up)
c       fact     :  EXP( - utaupr / ch ), where ch is the Chapman factor
c       plsorc   :  Planck source function (thermal)
c+---------------------------------------------------------------------+
      INTEGER layru(*)
      LOGICAL lyrcut, planck, prnt(*), spher
      REAL ch(*), cmu, dfdt(*), flup(*), fldir(*), fldn(*), kk( * ), 
     $     ll( mxcmu,* ), oprim(*), rfldir(*), rfldn(* ), rr(*), 
     $     ssalb(*), taucpr( 0:* ), u0c( mxcmu,mxulv ), uavg(*),
     $     utau(*), utaupr(*), xp_0(*), xp_1(*), yb_0d(*), yb_0u(*),
     $     yb_1d(*), yb_1u(*), yp_0d(*), yp_0u(*), yp_1d(*), yp_1u(*),
     $     zb_a(*), zp_a(*)
c
      IF ( prnt(2) )  WRITE( *,1010 )
c
c Zero twostr output arrays
c
      CALL  tzeroit( u0c, mxulv*mxcmu )
      CALL  tzeroit( rfldir, maxulv )
      CALL  tzeroit( fldir,  mxulv )
      CALL  tzeroit( rfldn,  maxulv )
      CALL  tzeroit( fldn,   mxulv )
      CALL  tzeroit( flup,   maxulv )
      CALL  tzeroit( uavg,   maxulv )
      CALL  tzeroit( dfdt,   maxulv )
c
c Loop over user levels
c
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
c
      DO 102  lu = 1, ntau
c
         lyu = layru(lu)
         IF ( lyrcut .AND. lyu.GT.ncut ) THEN
c
c No radiation reaches this level
c
            fdntot = 0.0
            fnet   = 0.0
            plsorc = 0.0
            GO TO 90                    ! ** Done with this layer
         END IF
c
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
c
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
c
c Calculate fluxes and mean intensities
c
c Downward and upward fluxes from Eq. KST(9)
c
         fldn( lu )  = 2.0 * pi * cmu * u0c(1,lu)
         flup( lu )  = 2.0 * pi *  cmu * u0c( 2,lu )
         fdntot      = fldn( lu ) + fldir( lu )
         fnet        = fdntot - flup( lu )
         rfldn( lu ) = fdntot - rfldir( lu )
c
c Mean intensity from Eq. KST(10)
c
         uavg(lu)   = u0c( 1,lu ) + u0c( 2, lu )
         uavg( lu ) = ( 2.0 * pi * uavg(lu) + dirint ) / ( 4.*pi )
c
c Flux divergence from Eqs. KST(11-12)
c
         plsorc = (1./(1.-oprim(lyu)))*EXP(-zp_a(lyu)*utaupr(lu))*
     $                          (xp_0(lyu) + xp_1(lyu)* utaupr(lu))
         dfdt( lu ) =  (1.0-ssalb(lyu)) * 4.*pi* (  uavg(lu) - plsorc )
c
 90      IF( prnt(2) )  WRITE( *,1020 ) utau(lu), lyu, rfldir(lu),
     $                                 rfldn(lu), fdntot, flup(lu),
     $                                 fnet, uavg(lu), plsorc, dfdt(lu)
 102  CONTINUE
c
 1010 FORMAT( //, 21X,
     $ '<----------------------- Fluxes ----------------------->', /,
     $ '   optical  compu    downward    downward    downward     ',
     $ ' upward                    mean      Planck   d(net flux)', /,
     $ '     depth  layer      direct     diffuse       total     ',
     $ 'diffuse         net   intensity      source   / d(op dep)', / )
 1020 FORMAT( F10.4, I7, 1P,7E12.3, E14.3 )
c
      RETURN
      END
c
      SUBROUTINE hopsol( biggest, ch, chtau, cmu, fbeam, ggprim, kk,
     $     ncut,nlyr, oprim, pi, pkag, pkagc, planck, radius, rr,
     $     smallest,spher,taucpr, umu0, xb_0d, xb_0u, xb_1d, xb_1u,
     $     xp_0,xp_1, yb_0d, yb_0u, yb_1d, yb_1u, yp_0d, yp_0u,
     $     yp_1d, yp_1u, zb_a, zp_a )
c
c Calculates the homogenous and particular solutions to the
c radiative transfer equation in the two-stream approximation,
c for each layer in the medium.
c
c    I n p u t     v a r i a b l e s:
c
c       ch       :  Chapman correction factor
c       cmu      :  Abscissae for gauss quadrature over angle cosine
c       ncut     :  Number of computational layer where absorption
c                     optical depth exceeds -abscut-
c       oprim    :  Delta-m scaled single scattering albedo
c       pkag,c   :  Planck function in each layer
c       spher    :  spher = true => spherical geometry invoked
c       taucpr   :  Cumulative optical depth (delta-m-scaled)
c       (remainder are 'twostr' input variables)
c
c   O u t p u t     v a r i a b l e s:
c
c       kk       :  Eigenvalues 
c       rr       :  Eigenvectors at polar quadrature angles 
c       xp_0,    :  Thermal source coefficients x-sub-zero and
c        xp_1         x-sub-one in Eq.  KST(22)
c       yb_0d,u, :  Thermal source vectors, Eqs. KST(24-25)
c        yb_1d,u     
c       yp_0d,u, :  Beam source vectors, Eq. KST(24-25)
c        yp_1d,u     
c       zb_a     :  Beam source coefficient alfa in Eq. KST(22)
c       zp_a     :  Thermal source coefficient alfa in Eq. KST(22)
c
c
      LOGICAL planck, spher
      REAL ch(*), chtau(*), ggprim(*), kk(*), large, oprim(*),
     $      pkag(0:*), pkagc(*), rr(*), taucpr(0:*), xb_0d(*),
     $     xb_0u(*), xb_1d(*), xb_1u(*), xp_0(*), xp_1(*), yb_0d(*),
     $     yb_0u(*), yb_1d(*), yb_1u(*), yp_0d(*), yp_0u(*), yp_1d(*),
     $     yp_1u(*), zb_a(*), zp_a(*)
c
c The calculation of the particular solutions require some care, small
c and big have been set so that no problems should occurr on 32-bits
c machine running single precision
c
      big   = sqrt(biggest) / 1.E+10
      small = 1.E+30*smallest
c
c ===================  Begin loop on computational layers  =============
c
      DO 100  LC = 1, NCUT
         ls = 2*(lc-1) + 1
c
c Calculate eigenvalues -kk- and eigenvector -rr-, Eqs. KST(20-21) 
c
         beta   = 0.5 * ( 1.-3.*ggprim(lc)*cmu*cmu )
         fact1  = 1. - oprim(lc)
         fact2  = 1. - oprim(lc) + 2.*oprim(lc)*beta 
         kk(lc) = (1./cmu) * sqrt( fact1 * fact2 )
         rr(lc) = ( sqrt(fact2)-sqrt(fact1) ) /
     $                                    ( sqrt(fact2)+sqrt(fact1) )
c
         IF ( fbeam.GT.0.0 ) THEN
c
c Set coefficients in KST(22) for beam source
c
           q_1 = (fbeam/(4.*pi))*oprim(lc)*(1.-3.*ggprim(lc)*cmu*umu0)
           q_2 = (fbeam/(4.*pi))*oprim(lc)*(1.+3.*ggprim(lc)*cmu*umu0) 
c
           IF ( umu0 .GE. 0.0 ) THEN
              qq = q_2
           ELSE
              qq = q_1
           ENDIF
c
           IF ( spher ) THEN
             q0a = EXP(-chtau(ls-1) )
             q0 = q0a*qq
             IF ( q0 .LE. small) THEN
                q1a = 0.0
                q2a = 0.0
             ELSE
                q1a = EXP(-chtau(ls) )
                q2a = EXP(-chtau(ls+1) )
             ENDIF
           ELSE IF ( .NOT. spher ) THEN
             q0a = EXP(-taucpr(lc-1)/umu0)
             q0 = q0a*qq
             IF ( q0 .LE. small) THEN
                q1a = 0.0
                q2a = 0.0
             ELSE
                q1a = EXP(- ( taucpr(lc-1)+taucpr(lc) )/ (2.*umu0) )
                q2a = EXP(-taucpr(lc)/umu0)
             ENDIF
           ENDIF
           q1 = q1a*qq
           q2 = q2a*qq
c
c Calculate alpha coefficient 
c
           deltat = taucpr(lc) - taucpr(lc-1)
           zb_a(lc) = 1./ch(lc)
           large = LOG(biggest)-20. 
           IF( ABS(zb_a(lc)*taucpr(lc-1)) .GT. large .OR.
     $          ABS(zb_a(lc)*taucpr(lc)) .GT. large ) zb_a(lc) = 0.0
c
c Dither alpha if it is close to an eigenvalue
c
           denomb =  fact1 * fact2 - (zb_a(lc)*cmu)*(zb_a(lc)*cmu)
           IF ( denomb .LT. 1.E-03 ) THEN
              zb_a(lc) = 1.02*zb_a(lc)
           ENDIF
           q0 = q0a * q_1
           q2 = q2a * q_1
c
c Set constants in Eq. KST(22)
c
           IF ( deltat .LT. 1.E-07 ) THEN
              xb_1d(lc) = 0.0
           ELSE
              xb_1d(lc) = (1./deltat)*(q2*EXP(zb_a(lc)*taucpr(lc)) 
     $             -q0*EXP(zb_a(lc)*taucpr(lc-1)))
           ENDIF
           xb_0d(lc) = q0 * EXP(zb_a(lc)*taucpr(lc-1)) -
     $          xb_1d(lc)*taucpr(lc-1)
           q0 = q0a * q_2
           q2 = q2a * q_2
           IF ( deltat .LT. 1.E-07 ) THEN
              xb_1u(lc) = 0.0
           ELSE
              xb_1u(lc) = (1./deltat)*(q2*EXP(zb_a(lc)*taucpr(lc)) 
     $             -q0*EXP(zb_a(lc)*taucpr(lc-1)))
           ENDIF
           xb_0u(lc) = q0 * EXP(zb_a(lc)*taucpr(lc-1)) -
     $          xb_1u(lc)*taucpr(lc-1)
c     
c Calculate particular solutions for incident beam source in 
c pseudo-spherical geometry, Eqs. KST(24-25)
c
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
c     
         ENDIF
c
         IF ( planck  ) THEN
c
c Set coefficients in KST(22) for thermal source
c
c Calculate alpha coefficient 
c
            small = 1.E+20*smallest
            q0 = (1.-oprim(lc)) * pkag(lc-1)
            q1 = (1.-oprim(lc)) * pkagc(lc)
            q2 = (1.-oprim(lc)) * pkag(lc)
            deltat = taucpr(lc) - taucpr(lc-1)
c
c Case 1: source small at bottom layer
c
            IF ( (q2.LT.(q0*1.E-02) .OR. q2.LE.small )
     $           .AND. q1.GT.small .AND. q0.GT.small ) THEN
c
c alpha Eq. KS(50)
c
               zp_a(lc) = (2./deltat) * LOG( q0/q1 )
               IF ( zp_a(lc) .GT. big )    zp_a(lc) = big
               IF ( zp_a(lc)*taucpr(lc-1) .GE. ALOG(big) ) THEN
                  xp_0(lc) =  big
               ELSE
                  xp_0(lc) = q0
               ENDIF
               xp_1(lc) = 0.0
c     
c Case 2: Source small at center and bottom of layer
c
            ELSE IF ( (q2.LE.(q1*1.E-02) .OR. q2.LE.small ) .AND.
     $              ((q1.LE.(q0*1.E-02)) .OR. (q1.LE.small))
     $              .AND. (q0.GT.small) ) THEN
c     
               zp_a(lc)  =   big / taucpr(ncut)
               xp_0(lc) = q0
               xp_1(lc) = 0.0
c     
c     Case 3:All sources zero
c
            ELSE IF ( q2.LE.small .AND. q1.LE.small
     $              .AND. q0.LE.small) THEN
               zp_a(lc)  = 0.0
               xp_0(lc) = 0.0
               xp_1(lc) = 0.0
c     
c     Case 4: Sources same at center, bottom and top of layer
c     or layer optically very thin
c
            ELSE IF ( (ABS((q2-q0)/q2).LT.1.E-04) .AND.
     $              (ABS((q2-q1)/q2).LT.1.E-04)
     $              .OR. deltat.LT. 1.E-04           ) THEN
c     
               zp_a(lc)  = 0.0
               xp_0(lc) = q0
               xp_1(lc) = 0.0
c     **  Case 5: Normal case
            ELSE
               arg = (q1/q2)**2. - q0/q2
               IF ( arg .LT. 0.0 ) arg = 0.0
c     
c alpha Eq. (44). For source that has its maximum at the top of the
c layer, use negative solution
c
               sgn = 1.0
               IF ( pkag(lc-1) .GT. pkag(lc) ) sgn = -1.
               fact3 = LOG(q1/q2 + sgn*SQRT(arg) )
               IF ( ABS(fact3) .LE. 0.005 ) THEN ! Be careful with log of
                  q1 = 0.99 * q1 ! numbers close to one
                  fact3 = LOG(q1/q2 + sgn*SQRT(arg) )
               ENDIF 
               zp_a(lc) = (2./deltat) * fact3
               IF(ABS(zp_a(lc)*taucpr(lc)) .GT. 
     $              (LOG(biggest)-LOG(q0*100.) ) )    zp_a(lc) = 0.0 
c     
c Dither alpha if it is close to an eigenvalue
c
               denomp =  fact1 * fact2 - (zp_a(lc)*cmu)*(zp_a(lc)*cmu)
               IF ( denomp .LT. 1.E-03 ) THEN
                  zp_a(lc) = 1.01*zp_a(lc)
               ENDIF
c
c Set constants in Eqs. KST(22)
c
               IF ( deltat .LT. 1.E-07 ) THEN
                  xp_1(lc) = 0.0
               ELSE
                  xp_1(lc) = (1./deltat)*(q2*EXP(zp_a(lc)*taucpr(lc)) 
     $                 -q0*EXP(zp_a(lc)*taucpr(lc-1)))
               ENDIF
               xp_0(lc) = q0 * EXP(zp_a(lc)*taucpr(lc-1)) -
     $              xp_1(lc)*taucpr(lc-1)
            ENDIF
c     
c Calculate particular solutions Eqs. KST(24-25) for internal thermal source
c
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
c     
         END IF
c     
 100  CONTINUE
c
c ===================  End loop on computational layers  ===============
c
      RETURN
      END
c
      SUBROUTINE tprtinp( albedo, btemp, deltam, dtauc, fbeam, fisot,   
     $     flyr, gg, lyrcut, nlyr, planck, nstr, ntau, oprim,
     $     spher, ssalb, tauc, taucpr, temis, temper, ttemp, umu0,
     $     utau, wvnmlo, wvnmhi )
c
c Print values of input variables
c
      LOGICAL  deltam, lyrcut, planck, spher
      REAL     dtauc(*), flyr(*), gg(*), oprim(*), ssalb(*),
     $         tauc( 0:* ), taucpr( 0:* ), temper( 0:* ), utau(*)
c
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
c
      RETURN
c
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
c
      END
c
      SUBROUTINE settwo( albedo, biggest, bplank, btemp, ch, chtau,
     $     cmu, deltam, 
     $     dtauc, dtaucpr, expbea, fbeam, flyr, gg, ggprim,
     $     layru, lyrcut, mxcly, ncut, nlyr, nn, nstr, ntau, oprim,
     $     pkag, pkagc, planck, radius, smallest, spher, ssalb, tauc, 
     $     taucpr, temis, temper, tplank, ttemp, umu0, usrtau,
     $     utau, utaupr, wvnmlo, wvnmhi, zd)
c
c Perform miscellaneous setting-up operations
c
c Routines called:  errmsg, tzeroit
c
c Input :  All are 'twostr' input variables (see doc file)
c
c Output:  ntau,utau   If usrtau = false
c          bplank      Intensity emitted from bottom boundary
c          ch          The Chapman factor
c          cmu         Computational polar angle
c          expbea      Transmission of direct beam
c          flyr        Truncated fraction in delta-m method
c          layru       Computational layer in which -utau- falls
c          lyrcut      Flag as to whether radiation will be zeroed
c                      below layer -ncut-
c          ncut        Computational layer where absorption
c                      optical depth first exceeds -abscut-
c          nn          nstr / 2  =  1
c          nstr        No. of streams (=2)
c          oprim       Delta-m-scaled single-scatter albedo
c          pkag,c      Planck function in each layer
c          taucpr      Delta-m-scaled optical depth
c          tplank      Intensity emitted from top boundary
c          utaupr      Delta-m-scaled version of -utau-
c
c Internal variables
c          abscut      Absorption optical depth, medium is cut off below
c                      this depth
c          tempc       Temperature at center of layer, assumed
c                      to be average of layer boundary temperatures
c
      LOGICAL  deltam, lyrcut, planck, spher, usrtau, first
      INTEGER  layru(*)
      REAL ch(*), chtau(*), cmu, dtauc(*), dtaucpr(*),
     $     expbea(0:*), flyr(*), gg(*), ggprim(*), oprim(*), pkag(0:*),
     $     pkagc(*), ssalb(*), tauc(0:*), taucpr(0:*),
     $     temper(0:*), utau(*), utaupr(*), zd(0:*)
c
      DATA  abscut / 10. /, first / .TRUE. /
c
      IF ( first ) THEN
         first    = .FALSE.
         smallest = r1mach(1)
         biggest  = r1mach(2)
         pi       = 2. * ASIN( 1.0 )
         nstr     = 2
         nn       = nstr / 2
      ENDIF
c
      IF ( .NOT.usrtau ) THEN
c
c Set output levels at computational layer boundaries
c
         ntau = nlyr + 1
         DO 30  lc = 0, ntau-1
            utau(lc+1) = tauc(lc)
 30      CONTINUE
      END IF
c
c Apply delta-m scaling and move description of computational
c layers to local variables
c
      expbea( 0 ) = 1.0
      zenang      = ACOS(umu0) * 180. / pi
      IF( spher .AND. umu0 .LT. 0.0 ) expbea(0) =
     $          EXP(-chapmn(1,0.0,tauc,nlyr, zd,dtauc,zenang,radius) )
      CALL  tzeroit( taucpr(0), mxcly+1 )
      CALL  tzeroit( expbea(1), mxcly )
      CALL  tzeroit( flyr, mxcly )
      CALL  tzeroit( oprim, mxcly )
      abstau = 0.0
      DO  60  lc = 1, nlyr
         IF ( abstau.LT.abscut )  ncut = lc
         abstau = abstau + ( 1. - ssalb(lc) ) * dtauc(lc)
c
         IF ( .NOT.deltam )  THEN
            oprim(lc)  = ssalb(lc)
            taucpr(lc) = tauc(lc)
            f          = 0.0
            ggprim(lc) = gg(lc)
            dtaucpr(lc)= dtauc(lc)
         ELSE
c
c Do delta-m transformation Eqs. WW(20a,20b,14)
c
            f = gg(lc) * gg(lc)
            taucpr(lc) = taucpr(lc-1) + ( 1. - f*ssalb(lc) ) * dtauc(lc)
            oprim(lc)  = ssalb(lc) * ( 1. - f ) / ( 1. - f * ssalb(lc) )
            ggprim(lc) =  (gg(lc)-f) / (1.-f)
            dtaucpr(lc)= taucpr(lc) - taucpr(lc-1)
         ENDIF
c
         flyr(lc)   = f
c
 60   CONTINUE
c
c If no thermal emission, cut off medium below absorption optical
c depth = abscut ( note that delta-m transformation leaves absorption
c optical depth invariant ).  Not worth the trouble for one-layer
c problems, though.
c
      lyrcut = .FALSE.
      IF ( abstau.GE.abscut .AND. .NOT. planck
     $                           .AND. nlyr.GT.1 )  lyrcut =.TRUE.
      IF ( .NOT.lyrcut )  ncut = nlyr
c
c Calculate chapman function is spherical geometry, set expbea and ch
c for beam source.
c
      IF ( fbeam.GT.0.0 )  THEN
         chtau(0) = 0.0
         DO lc = 1, ncut
            expbea(lc) = 0.0
            IF ( spher ) THEN
               ls = 2*(lc-1) + 1
               taup   = taucpr(lc-1) + dtaucpr(lc)/2.0
               chtau(ls) = chapmn(lc,taup,taucpr,nlyr,
     $              zd,dtaucpr,zenang,radius)
               chtau(ls+1) = chapmn(lc,taucpr(lc),taucpr,nlyr,
     $              zd,dtaucpr,zenang,radius)
               ch(lc) = taup/chtau(ls)
               expbea(lc) = EXP(-chtau(ls+1) )
            ELSE IF ( .NOT. spher ) THEN
               ch(lc)     = umu0
               expbea(lC) = EXP( - taucpr(lc) / umu0 )
            ENDIF
         ENDDO
      ENDIF
c
c Set arrays defining location of user output levels within 
c delta-m-scaled computational mesh
c 
      DO 90  lu = 1, ntau
         DO 70 lc = 1, nlyr
            IF ( utau(lu).GE.tauc(lc-1) .AND. utau(lu).LE.tauc(lc) )
     $           GO TO 80
 70      CONTINUE
         lc = nlyr
c
 80      utaupr(lu) = utau(lu)
         IF(deltam) utaupr(lu) = taucpr(lc-1) + (1.-ssalb(lc)*flyr(lc))
     $                                        * (utau(lu) - tauc(lc-1))
         layru(lu) = lc
 90   CONTINUE
c
c Set computational polar angle cosine for double gaussian 
c quadrature; cmu = 0.5, or  single gaussian quadrature; cmu = 1./sqrt(3)
c See KST for discussion of which is better for your specific application
c
      IF ( planck .AND. fbeam .EQ. 0.0 ) THEN
         cmu =  0.5
      ELSE
         cmu = 1./SQRT(3.0)
      ENDIF
c
c Calculate planck functions
c
      IF ( .NOT. planck )  THEN
         bplank = 0.0
         tplank = 0.0
         CALL  tzeroit( pkag, mxcly+1 )
         CALL  tzeroit( pkagc, mxcly )
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
c
      SUBROUTINE solbcd( albedo, b, bplank, cband, cmu, diag, expbea,
     $     fbeam, fisot, ipvt, kk, ll, lyrcut, mi, mi9m2, mxcmu,
     $     ncol, ncut, nlyr, nn, nnlyri, nstr, pi, rr, subd,
     $     superd, sol,tplank, taucpr, umu0, yb_0d, yb_0u, yb_1d,
     $     yb_1u, yp_0d, yp_0u, yp_1d, yp_1u, zb_a, zp_a )
c
c Construct right-hand side vector -b- for general boundary conditions
c and solve system of equations obtained from the boundary conditions
c and the continuity-of-intensity-at-layer-interface equations.
c
c Routines called: sgbfa, sgbsl, tzeroit
c
c I n p u t      v a r i a b l e s:
c
c       bplank   :  Bottom boundary thermal emission
c       cband    :  Left-hand side matrix of linear system Eqs. KST(38-41)
c                   in banded form required by linpack solution routines
c       cmu      :  Abscissae for gauss quadrature over angle cosine
c       expbea   :  Transmission of incident beam, EXP(-taucpr/ch)
c       lyrcut   :  Logical flag for truncation of comput. layer
c       ncol     :  Counts of columns in -cband-
c       nn       :  Order of double-gauss quadrature (nstr/2)
c       ncut     :  Total number of computational layers considered
c       nstr     :  No. of streams (=2)
c       tplank   :  Top boundary thermal emission
c       taucpr   :  Cumulative optical depth (delta-m-scaled)
c       yb_0d,u, :  Thermal source vectors, Eq. KST(24-25)
c        yb_1d,u     
c       yp_0d,u, :  Beam source vectors, Eq. KST(24-25)
c        yp_1d,u     
c       zb_a     :  Beam source coefficient alfa in Eq. KST(22)
c       zp_a     :  Thermal source coefficient alfa in Eq. KST(22)
c       (Remainder are 'twostr' input variables)
c
c O u t p u t     v a r i a b l e s:
c
c       b        :  Right-hand side vector of Eqs. KST(38-41) going into
c                   *sgbsl*; returns as solution vector of Eqs. KST(38-41)
c                   constants of integration 
c      ll        :  Permanent storage for -b-, but re-ordered
c
c I n t e r n a l    v a r i a b l e s:
c
c       it       :  Pointer for position in  -b-
c       ncd      :  Number of diagonals below or above main diagonal
c       rcond    :  Indicator of singularity for -cband-
c       z        :  Scratch array required by *sgbco*
c+---------------------------------------------------------------------+
c
      INTEGER  ipvt(*)
      LOGICAL  lyrcut
      REAL b(*), cband( mi9m2,nnlyri ), cmu, diag(*), expbea(0:*),
     $     kk(*), ll( mxcmu,* ), rr(*), subd(*), superd(*), sol(*),
     $     taucpr( 0:* ), yb_0d(*), yb_0u(*), yb_1d(*), yb_1u(*),
     $     yp_0d(*), yp_0u(*), yp_1d(*), yp_1u(*), zb_a(*), zp_a(*)
c
c First top row, top boundary condition
c
      irow        = 1
      lc          = 1
c     subd(irow)  = undefined
      diag(irow)  = rr(lc) * EXP(-kk(lc) * taucpr(lc))
      superd(irow)= 1.0
c
c next from layer no. 2 to nlyr -1
c
      nloop = ncut - 1 
      DO lc = 1, nloop
         irow         = irow + 1
         wk0          = EXP(-kk(lc) * (taucpr(lc) - taucpr(lc-1)))
         wk1          = EXP(-kk(lc+1) * (taucpr(lc+1) - taucpr(lc)))
         subd(irow)   = 1.0 - rr(lc) * rr(lc+1)
         diag(irow)   = ( rr(lc) -  rr(lc+1 ) ) * wk0
         superd(irow) = - ( 1. - rr(lc+1) * rr(lc+1 ) ) * wk1
         irow         = irow + 1
         subd(irow)   = ( 1.0 - rr(lc) * rr(lc) ) * wk0 
         diag(irow)   = ( rr(lc) -  rr(lc+1 ) ) * wk1
         superd(irow) = - ( 1. - rr(lc+1) * rr(lc ) ) 
      ENDDO
c
c bottom layer
c
      irow         = irow + 1
      lc           = ncut
c     superd(irow) = undefined
      wk           = EXP( -kk(lc) * (taucpr(lc) - taucpr(lc-1)) )
      IF ( lyrcut ) THEN
         subd(irow) = 1.0
         diag(irow) = rr(lc) * wk
      ELSE
         subd(irow) = 1.0 - 2.0*albedo*cmu*rr(lc)
         diag(irow) = ( rr(lc) - 2.0*albedo*cmu ) * wk         
      ENDIF
c
      CALL  tzeroit( b, nnlyri )
c
c Construct -b-,  for parallel beam + bottom reflection +
c thermal emission at top and/or bottom
c
c Top boundary, right-hand-side of Eq. KST(28)
c
      lc      = 1
      irow    = 1
      b(irow) = - yb_0d(lc) - yp_0d(lc) +fisot +tplank
c
c Continuity condition for layer interfaces,
c right-hand-side of Eq. KST(29)
c
      DO   lc = 1, nloop
         rpp1_m = EXP(-zb_a(lc+1)*taucpr(lc))*
     $        (yb_0d(lc+1)+yb_1d(lc+1)*taucpr(lc)) 
     $        + EXP(-zp_a(lc+1)*taucpr(lc))*
     $        (yp_0d(lc+1)+yp_1d(lc+1)*taucpr(lc))
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
      ENDDO
c
c Bottom boundary
c
      irow = irow + 1
      lc   = ncut
      IF ( lyrcut ) THEN
c
c Right-hand-side of Eq. KST(30)
c
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
c
c solve for constants of integration by inverting matrix KST(38-41)
c
      nrow = irow
c
          CALL tzeroit( cband, mi9m2*nnlyri)
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
c
          CALL sgbfa(cband, mi9m2, nrow, 1, 1, ipvt, info )
          job = 0
          CALL sgbsl(cband, mi9m2, nrow, 1, 1, ipvt, b, job )
c
c unpack
c
          irow = 1
          DO lc = 1, ncut
             ll(1,lc) = b(irow)      ! downward direction
             irow     = irow + 1
             ll(2,lc) = b(irow)      ! upward direction
             irow     = irow + 1
          ENDDO
c
      RETURN
      END
c
      SUBROUTINE tzeroal( ch, kk, ll, mxcmu, mxcly,
     $     nnlyri, rr, xb_0d, xb_0u, xb_1d, xb_1u,
     $     xp_0, xp_1, yb_0d, yb_0u, yb_1d, yb_1u, yp_0d,
     $     yp_0u, yp_1d, yp_1u, zb_a, zp_a )
c
c Zero arrays
c
      REAL  ch(*), kk(*), ll(mxcmu,*), rr(*), 
     $      xb_0d(*), xb_0u(*), xb_1d(*),
     $      xb_1u(*), xp_0(*), xp_1(*), yb_0d(*), yb_0u(*),
     $      yb_1d(*), yb_1u(*), yp_0d(*), yp_0u(*), yp_1d(*), 
     $      yp_1u(*), zb_a(*), zp_a(*)
c
      CALL  tzeroit( ch    , mxcly )
      CALL  tzeroit( kk,     mxcly )
      CALL  tzeroit( ll,     mxcmu*mxcly )
      CALL  tzeroit( rr,     mxcly )
      CALL  tzeroit( xb_0d, mxcly )
      CALL  tzeroit( xb_0u, mxcly )
      CALL  tzeroit( xb_1d, mxcly )
      CALL  tzeroit( xb_1u, mxcly )
      CALL  tzeroit( xp_0, mxcly )
      CALL  tzeroit( xp_1, mxcly )
      CALL  tzeroit( yb_0d, mxcly )
      CALL  tzeroit( yb_0u, mxcly )
      CALL  tzeroit( yb_1d, mxcly )
      CALL  tzeroit( yb_1u, mxcly )
      CALL  tzeroit( yp_0d, mxcly )
      CALL  tzeroit( yp_0u, mxcly )
      CALL  tzeroit( yp_1d, mxcly )
      CALL  tzeroit( yp_1u, mxcly )
      CALL  tzeroit( zb_a, mxcly )
      CALL  tzeroit( zp_a, mxcly )
c
      RETURN
      END
c
      SUBROUTINE tzeroit( A, LENGTH )
C
C         ZEROS A REAL ARRAY -A- HAVING -LENGTH- ELEMENTS
C
      REAL  A(*)
C
      DO 10  L = 1, LENGTH
         A( L ) = 0.0
10    CONTINUE
C
      RETURN
      END
      REAL FUNCTION  TPLKAVG ( WNUMLO, WNUMHI, T )
C
C        COMPUTES PLANCK FUNCTION INTEGRATED BETWEEN TWO WAVENUMBERS,
c        except if wnmulo .EQ. wnmuhi, then the Planck function at 
c        wnumlo is returned
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
      INTEGER  SMALLV
      REAL     C2, CONC, D(2), EPSIL, EX, MV, P(2), SIGMA, SIGDPI,
     $         V(2), VCUT, VCP(7), VSQ
      SAVE     CONC, VMAX, EPSIL, SIGDPI
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
c
      IF ( wnumhi .eq. wnmulo ) THEN
         wvn  =  wnumhi
         arg  = EXP( - C2 * wvn / T )
         plkavg = c1 * (wvn**3.) * arg / ( 1. - arg )
         RETURN
      ENDIF
c     
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
