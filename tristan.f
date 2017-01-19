      subroutine trist(tin,pin,hminus,h2minus)
c     Returns ratio of H-/H2 number densities as a function of t and p
      implicit double precision (a-h,o-z)
      integer pn,pe
      parameter (pn=51,pe=20)
      common /tristan/ p(pn),t(pn),ro(pn,pn),xnel(pn,pn),ab(pe,pn,pn),
     &     code(pe),np,nt
c Here is a correspondance for the codes:
c  1.			H
c  1.01			H+
c  2.			He
c           100.00	H-
c           101.00	H2
c           608.00	CO
c           707.00	N2
c           822.00	TiO
c         10108.00	H2O
c         10116.00	H2S
c         10607.00	HCN
c       1010107.00	NH3
c     101010106.00	CH4

	   tl = dlog10(tin)
 	   pl = dlog10(pin*1e6)

      CALL LOCATE (t,nt,tl,kt)  
      CALL LOCATE (p,np,pl,kp)  
      FACTkt= (-T(Kt)+Tl)/(T(Kt+1)-T(Kt))
      FACTkp= (-P(Kp)+Pl)/(P(Kp+1)-P(Kp))
 	  if (kp.eq.0) then
c       we are at low pressure, use the lowest pressure point
 	    factkp = 0.0
 		kp = 1
   	  endif
 	  if (kp.eq.np) then
c       we are at high pressure, use the highest pressure point
 	    factkp = 0.0
   	  endif
 	  if (kt.eq.0) then
c       we are at low temp, use the lowest temp point
 	    factkt = 0.0
 		kt = 1
   	  endif
 	  if (kt.eq.nt) then
c       we are at high temp, use the highest temp point
 	    factkt = 0.0
   	  endif
	y1 = ab(4,kt,kp)-ab(5,kt,kp)	
        y2 = ab(4,kt+1,kp)-ab(5,kt+1,kp)
        y3 = ab(4,kt+1,kp+1)-ab(5,kt+1,kp+1)
        y4 = ab(4,kt,kp+1)-ab(5,kt,kp+1)

c	hminus = [H-]/[H2]
	hminus = (1.d0-factkt)*(1.d0-factkp)*y1 + 
     1   factkt*(1.d0-factkp)*y2 + factkt*factkp*y3 + (1.d0-factkt)*factkp*y4
	hminus = 10.d0**hminus

	y1 = xnel(kt,kp)-ab(5,kt,kp)	
        y2 = xnel(kt+1,kp)-ab(5,kt+1,kp)
        y3 = xnel(kt+1,kp+1)-ab(5,kt+1,kp+1)
        y4 = xnel(kt,kp+1)-ab(5,kt,kp+1)

c	h2minus = [electrons]/[H2]
	h2minus = (1.d0-factkt)*(1.d0-factkp)*y1 + 
     1   factkt*(1.d0-factkp)*y2 + factkt*factkp*y3 + (1.d0-factkt)*factkp*y4
	h2minus = 10.d0**h2minus

      RETURN
      END


      subroutine trist1(tin,pin,ch4x,xnh3x,xh2o,xn2x,xco)
c     Returns ratio of []/H2 number densities as a function of t and p
      implicit double precision (a-h,o-z)
      integer pn,pe
      parameter (pn=51,pe=20)
      common /tristan/ p(pn),t(pn),ro(pn,pn),xnel(pn,pn),ab(pe,pn,pn),
     &     code(pe),np,nt
      COMMON /metals/ f_metal,xnh3save,xh2save,xh2osave,xtiosave,xfehsave
c Here is a correspondance for the codes:
c  1.			H
c  1.01			H+
c  2.			He
c           100.00	H-
c           101.00	H2
c           608.00	CO
c           707.00	N2
c           822.00	TiO
c         10108.00	H2O
c         10116.00	H2S
c         10607.00	HCN
c       1010107.00	NH3
c     101010106.00	CH4


	   tl = dlog10(tin)
 	   pl = dlog10(pin*1e6)

      CALL LOCATE (t,nt,tl,kt)  
      CALL LOCATE (p,np,pl,kp)  
      FACTkt= (-T(Kt)+Tl)/(T(Kt+1)-T(Kt))
      FACTkp= (-P(Kp)+Pl)/(P(Kp+1)-P(Kp))
 	  if (kp.eq.0) then
c       we are at low pressure, use the lowest pressure point
 	    factkp = 0.0
 		kp = 1
   	  endif
 	  if (kp.eq.np) then
c       we are at high pressure, use the highest pressure point
 	    factkp = 0.0
   	  endif
 	  if (kt.eq.0) then
c       we are at low temp, use the lowest temp point
 	    factkt = 0.0
 		kt = 1
   	  endif
 	  if (kt.eq.nt) then
c       we are at high temp, use the highest temp point
 	    factkt = 0.0
   	  endif

c	xch4 = [CH4]/[H2]
	y1 = ab(13,kt,kp)-ab(5,kt,kp)	
        y2 = ab(13,kt+1,kp)-ab(5,kt+1,kp)
        y3 = ab(13,kt+1,kp+1)-ab(5,kt+1,kp+1)
        y4 = ab(13,kt,kp+1)-ab(5,kt,kp+1)

	ch4x = (1.d0-factkt)*(1.d0-factkp)*y1 + 
     1   factkt*(1.d0-factkp)*y2 + factkt*factkp*y3 + (1.d0-factkt)*factkp*y4
	ch4x = 10.d0**ch4x

c	xnh3 = [NH3]/[H2]
	y1 = ab(12,kt,kp)-ab(5,kt,kp)	
        y2 = ab(12,kt+1,kp)-ab(5,kt+1,kp)
        y3 = ab(12,kt+1,kp+1)-ab(5,kt+1,kp+1)
        y4 = ab(12,kt,kp+1)-ab(5,kt,kp+1)

	xnh3 = (1.d0-factkt)*(1.d0-factkp)*y1 + 
     1   factkt*(1.d0-factkp)*y2 + factkt*factkp*y3 + (1.d0-factkt)*factkp*y4
	xnh3x = 10.d0**xnh3

c	xh2o = [H2O]/[H2]
	y1 = ab(9,kt,kp)-ab(5,kt,kp)	
        y2 = ab(9,kt+1,kp)-ab(5,kt+1,kp)
        y3 = ab(9,kt+1,kp+1)-ab(5,kt+1,kp+1)
        y4 = ab(9,kt,kp+1)-ab(5,kt,kp+1)

	xh2o = (1.d0-factkt)*(1.d0-factkp)*y1 + 
     1   factkt*(1.d0-factkp)*y2 + factkt*factkp*y3 + (1.d0-factkt)*factkp*y4
	xh2o = 10.d0**xh2o

c	xn2x = [N2]/[H2]
	y1 = ab(7,kt,kp)-ab(5,kt,kp)	
        y2 = ab(7,kt+1,kp)-ab(5,kt+1,kp)
        y3 = ab(7,kt+1,kp+1)-ab(5,kt+1,kp+1)
        y4 = ab(7,kt,kp+1)-ab(5,kt,kp+1)

	xn2x = (1.d0-factkt)*(1.d0-factkp)*y1 + 
     1   factkt*(1.d0-factkp)*y2 + factkt*factkp*y3 + (1.d0-factkt)*factkp*y4
	xn2x = 10.d0**xn2x

c	xco  = [CO]/[H2]
	y1 = ab(6,kt,kp)-ab(5,kt,kp)	
        y2 = ab(6,kt+1,kp)-ab(5,kt+1,kp)
        y3 = ab(6,kt+1,kp+1)-ab(5,kt+1,kp+1)
        y4 = ab(6,kt,kp+1)-ab(5,kt,kp+1)

	xco  = (1.d0-factkt)*(1.d0-factkp)*y1 + 
     1   factkt*(1.d0-factkp)*y2 + factkt*factkp*y3 + (1.d0-factkt)*factkp*y4
	xco  = 10.d0**xco 

      ch4x = f_metal * ch4x
	  xnh3x = f_metal * xnh3x
	  xh2o = f_metal * xh2o
	  xn2x = f_metal * xn2x
	  xco = f_metal * xco

      RETURN
      END




c  Here is the subroutine I use to calculate the bound-free absorption of
c  the hydrogen ion. The reference is in my paper "Are the giant planets..."
c  The result is given in cm2. In order to calculate an opacity in cm-1,
c  multiply by the number of H- ions/cm3. To have it in cm2/g, divide by
c  the density in g/cm3.
c  Caution: sbf_hm does not include stimulated emission.

      subroutine opa_hmbf(sigma,sbf_hm)
	implicit double precision (a-h,o-z)
c     Version: 05/10/92
c     Section efficace d'absorption bound-free de H- d'apres John (1988)
c     Auteur: T.Guillot
c Entree:
c     sigma: nombre d'onde [cm-1]
c Sortie:
c     sbf_hm: coefficient d'absorption [cm2]
	dimension coeff(6)
      double precision lambda,lambda0
      data coeff/152.519d0,49.534d0,-118.858d0,92.536d0,
     &     -34.194d0,4.982d0/
      data lambda0/1.6419d0/

      if (sigma.gt.1d4/lambda0) then
       lambda=1d4/sigma
       x=sqrt(1/lambda-1/lambda0)
       f=0.d0
       do i=6,1,-1
        f=f*x+coeff(i)
       enddo
       sbf_hm=(lambda*x)**3*f*1d-18
      else
       sbf_hm=0.d0
      endif

      return
      end



c**********************************************************************
c  Opacity for h2 minus
c  The result is given in cm4/dyn, so you should multiply that by
c	nh2*ne*k*T
c to obtain an opacity in cm-1.
c (nh2 is the number of h2 molecules/cm3, ne is the same for electrons,
c k is the boltzmann constant). 
c Unlike the other subroutine I sent you, the result *includes* the
c contribution due to stimulated emission.
c     Version: 30/12/95


      subroutine opa_tab_h2mff(t,sigma,sff_h2m)
c     Version: 13/04/93
c     Section efficace d'absorption free-free de H2- d'apres Bell (1980)
c     Auteur: T.Guillot
c Entrees:
c     t: temperature [K]
c     sigma: nombre d'onde [cm-1]
c Sortie:
c     sff_h2m: coefficient d'absorption [cm4/dyn]
      implicit none
      save
      double precision t,sigma,sff_h2m

      integer pnx,pny,spx,spy,pxt,pyt,pa,dims
      parameter (pnx=18, pny=9, spx=4, spy=4)
      parameter (pxt=2*spx+pnx, pyt=2*spy+pny, pa=pnx*pny, 
     &     dims=pnx*pny*spx*spy)

      integer nx,my,knotx,nxr,i,ny,mx,knoty,nyr,lx,ly

      double precision x(pnx),xt(pxt),xr(pnx),f(pa),s(dims),
     &     y(pny),yt(pyt),yr(pny)
      double precision xl,yl,fl,dfdx,dfdy

      logical init
      data init/.true./

      data (y(i),i=1,9)
     &     /0.5d0,0.8d0,1.0d0,1.2d0,1.6d0,2.0d0,2.8d0,3.6d0,10.d0/

c Les "-" sont dus au fait que x(i) doit etre une suite croissante
      data (x(i),i=1,18)
     &     /-151883d0,-113913d0,-91130d0,-60753d0,-45565d0,-36452d0,
     &     -30377d0,-22783d0,-18226d0,-15188d0,-11391d0,-9113d0,
     &     -7594d0,-6509d0,-5696d0,-5063d0,-4142d0,-3505d0/

      data (f(i),i=1,18)
     &     /7.16d1,4.03d1,2.58d1,1.15d1,6.47d0,4.15d0,2.89d0,1.63d0,
     &      1.05d0,7.36d-1,4.20d-1,2.73d-1,1.92d-1,1.43d-1,1.10d-1,
     &      8.70d-2,5.84d-2,4.17d-2/
      data (f(i),i=19,36)
     &     /9.23d1,5.20d1,3.33d1,1.48d1,8.37d0,5.38d0,3.76d0,2.14d0,
     &      1.39d0,9.75d-1,5.64d-1,3.71d-1,2.64d-1,1.98d-1,1.54d-1,
     &      1.24d-1,8.43d-2,6.10d-2/
      data (f(i),i=37,54)
     &     /1.01d2,5.70d1,3.65d1,1.63d1,9.20d0,5.92d0,4.14d0,2.36d0,
     &      1.54d0,1.09d0,6.35d-1,4.22d-1,3.03d-1,2.30d-1,1.80d-1,
     &      1.46d-1,1.01d-1,7.34d-2/
      data (f(i),i=55,72)
     &     /1.08d2,6.08d1,3.90d1,1.74d1,9.84d0,6.35d0,4.44d0,2.55d0,
     &      1.66d0,1.18d0,6.97d-1,4.67d-1,3.39d-1,2.59d-1,2.06d-1,
     &      1.67d-1,1.17d-1,8.59d-2/
      data (f(i),i=73,90)
     &     /1.18d2,6.65d1,4.27d1,1.91d1,1.08d1,6.99d0,4.91d0,2.84d0,
     &      1.87d0,1.34d0,8.06d-1,5.52d-1,4.08d-1,3.17d-1,2.55d-1,
     &      2.10d-1,1.49d-1,1.11d-1/
      data (f(i),i=91,108)
     &     /1.26d2,7.08d1,4.54d1,2.04d1,1.16d1,7.50d0,5.28d0,3.07d0,
     &      2.04d0,1.48d0,9.09d-1,6.33d-1,4.76d-1,3.75d-1,3.05d-1,
     &      2.53d-1,1.82d-1,1.37d-1/
      data (f(i),i=109,126)
     &     /1.38d2,7.76d1,4.98d1,2.24d1,1.28d1,8.32d0,5.90d0,3.49d0,
     &      2.36d0,1.74d0,1.11d0,7.97d-1,6.13d-1,4.92d-1,4.06d-1,
     &      3.39d-1,2.49d-1,1.87d-1/
      data (f(i),i=127,144)
     &     /1.47d2,8.30d1,5.33d1,2.40d1,1.38d1,9.02d0,6.44d0,3.90d0,
     &      2.68d0,2.01d0,1.32d0,9.63d-1,7.51d-1,6.09d-1,5.07d-1,
     &      4.27d-1,3.16d-1,2.40d-1/
c !!! Extrapolation lineaire pour T<1400K
      data (f(i),i=145,162)
     &     /2.19d2,1.26d2,8.13d1,3.68d1,2.18d1,1.46d1,1.08d1,7.18d0,
     &      5.24d0,4.17d0,3.00d0,2.29d0,1.86d0,1.55d0,1.32d0,1.13d0,
     &      8.52d-1,6.64d-1/

 2000 format((1x,1p8d10.3))

      if (init) then
       init=.false.
c       write(*,100)
c 100   format(1x,'Table des absorptions H2- free-free d''apres Bell')
       
       nx=pnx                   !nb de x
       ny=pny                   !nb de y
       mx=4                     !ordre des splines
       my=4
	
       call pp2s(x,xt,xr,nx,nxr,mx,knotx,f,s,
     &      y,yt,yr,ny,nyr,my,knoty)
      endif
       
c     interpolation par b-spline 2-D d'ordre (mx,my) au point (rol,tl)
      
      if (sigma.eq.0.) then
       xl=x(1)
      else
       xl=-1d8/sigma
      endif
      yl=5040.4/t
      
      if (xl.lt.x(1)) then
       call pp2d(x,xt,xr,nx,nxr,mx,x(1),knotx,lx,dfdx,f,s,fl,
     &      .true.,y,yt,yr,ny,nyr,my,yl,knoty,ly,dfdy)
       sff_h2m=fl*(xl/x(1))**2*1d-26	!d'apres Gould (1985)
      else if (xl.gt.x(nx)) then
       call pp2d(x,xt,xr,nx,nxr,mx,x(nx),knotx,lx,dfdx,f,s,fl,
     &      .true.,y,yt,yr,ny,nyr,my,yl,knoty,ly,dfdy)
       sff_h2m=fl*(xl/x(nx))*1d-26	!d'apres Gould (1985)
      else
       call pp2d(x,xt,xr,nx,nxr,mx,xl,knotx,lx,dfdx,f,s,fl,
     &      .true.,y,yt,yr,ny,nyr,my,yl,knoty,ly,dfdy)
       sff_h2m=fl*1d-26
      endif
      
c     rol ou tl sort de la table initiale ?
      
      if(yl .lt. y(1))then
       write(*,225)yl,y(1)
 225   format(1x,'!Sortie de table dans OPA_TAB_H2MFF: y=',1pd14.7,
     &      ' < ymin=',1pd14.7)
      elseif(yl .gt. y(ny))then
       write(*,226)yl,y(ny)
 226   format(1x,'!Sortie de table dans OPA_TAB_H2MFF: y=',1pd14.7,
     &      ' > ymax=',1pd14.7)
      endif

      return
      end
      
c************************************************************************
      
      subroutine pp2s(x,xt,xr,nx,nrx,mx,knotx,f,s,
     &     y,yt,yr,ny,nry,my,knoty)
c     Version: 01/08/95
      
c     Calcul des coefficients des PP pour interpolation 2D
c     Correspond a la premiere partie de PP2D et de PP2DD
c     Evite de reserver la memoire pour un tres grand tableau des ax
      
c     interpolation dans le tableau f(nx,ny) par produit tensoriel
c     de splines polynomiales d'ordres mx{y}>1, au point
c     (xx,yy) : x(1).le.xx.le.x(nx), y(1).le.yy.le.y(ny), on aura aussi
c      xt(lx).le.xx.lt.xt(lx+1)       yt(ly).le.yy.lt.yt(lt+1)
c     en entree prendre lx{y}=mx{y}
      
c     calul de knotx=2mx+nx, formation du tableau xt(knotx). idem pour y
c     !! ATTENTION le tableau f(nx,ny) initial est remplace par le
c     tableau des coefficients des B-splines en xy!!
c     ce calcul etant fait selon de Boor p.342
      
c     ax et ay doivent etre dimensionnes au moins a   max(nx*nx,ny*ny)
c     l et m  "       "       "       "       "       max(nx,ny)
      
c     les tables de travail xr, yr sont conservees car on peut faire appel
c     sequentiellement a des tables diverses
      
c     Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c     Version: 04/01/95
      
      implicit none
      
      integer dimab,dimxy,dimaxy,dimqxy,dimdxy
      parameter (dimab=501, dimxy=4)
      parameter (dimaxy=dimab**2, dimqxy=dimxy+1, dimdxy=dimxy**2)
      
      integer nx,ny,knotx,knoty,lx,ly,mx,my,l(dimab),m(dimab),i,j,k,
     &     llx,lly,ix,iy,sx,sy,nrx,nry
      
      double precision x(*),y(*),f(*),xt(*),yt(*),s(*),xr(*),yr(*),
     &     ax(dimaxy),ay(dimaxy),qx(dimqxy),qy(dimqxy),dx(dimdxy),
     &     dy(dimdxy),pivot,deriv,fac(dimxy)
      
      
 2000 format(1x,1p8e10.3)
      
      if(min(mx,my) .lt. 2)then
       write(6,*)'avec PP2$ il faut min(mx,my) > 2 on a: mx,my=',mx,my
       stop
      elseif(max(mx,my) .gt. dimxy)then
       write(6,*)'avec PP2$ il faut max(mx,my)<dimxy on a: mx,my,',
     &      'dimxy=',mx,my,dimxy
       write(6,*)'modifier le parametre dimxy'
       stop
      elseif(max(nx,ny) .gt. dimab)then
       write(6,*)'dans PP2$ max(nx,ny) est trop grand nx,ny',nx,ny
       write(6,*)'modifier le parametre dimab'
       stop
      endif
      
      
      do i=1,dimaxy
       ax(i)=0.
       ay(i)=0.
      end do
      
      lx=mx
      ly=my
      
c     formation des tableaux xt, yt et points de raccord
      
      call snoein(x,xt,nx,mx,knotx)
      call snoein(y,yt,ny,my,knoty)
      nrx=nx-mx+2
      nry=ny-my+2
      do ix=1,nrx
       xr(ix)=xt(ix+mx-1)
      enddo
      do iy=1,nry
       yr(iy)=yt(iy+my-1)
      enddo
      
c     ay <- matrice 1-D des interpolations en y
      
      do j=1,ny
       call slinf(y(j),yt,knoty,ly)
       call bval(y(j),yt,my,ly,qy)
       do k=1,my
        ay(ny*(ly-my+k-1)+j)=qy(k)
       enddo                    !k
      enddo                     !j
      
c     ax <- transposee de ay, inversion de ax et ay=f*ax
      
      call strans(ay,ny,ax)
      
      call minv(ax,ny,pivot,l,m)
      
      if(pivot .eq. 0.)then
c      write(6,*)'dans l''inversion 1 de ax'
       stop
      else
       call smatpr(f,nx,ny,ax,ny,ay)
      endif
      
c     ax <- matrice 1-D des interpolations en x
      
      do i=1,dimaxy
       ax(i)=0.
      enddo
      do j=1,nx
       call slinf(x(j),xt,knotx,lx)
       call bval(x(j),xt,mx,lx,qx)
       do k=1,mx
        ax(nx*(lx-mx+k-1)+j)=qx(k)
       enddo                    !k
      enddo                     !j
      
c     ax=inverse de ax et f=ax*ay
      
      call minv(ax,nx,pivot,l,m)
      if(pivot .eq. 0.)then
       write(6,*)'dans l''inversion 2 de ax'
       stop
      else
       call smatpr(ax,nx,nx,ay,ny,f)
      endif
      
c     calcul des coefficients des pp. : 12-21 de schumaker
      
c     de facon a optimiser l'algoritme de horner avec derivee
c     dans s, les coefficients du polynome (derives/fac(.)) sont classes
c     par puissance decroissante :
c     fy3x3, fy3x2, fy3x, fy3; fy2x3, fy2x2, fy2x, fy2;
c     fyx3, fyx2, fyx,fy; fx3, fx2, fx, f
      
c     s(ix,iy,lx,ly)=s(mx*(my*((nrx-1)*(ly-1)+lx-1)+iy-1)+ix)
c     il faut ajouter 1 a l'ordre de la derivation, ainsi
c     pour fy2x on aura ix=mx-2, iy=my-3 et, pour l'ordre de derivation i,j
c     correspondant aux indices k=i+1, l=j+1 pour bvald, on aura
c     ix=mx-k+1, iy=my-l+1
      
      fac(1)=1.                 !calcul de (i-1)!
      fac(2)=1.
      do i=3,dimxy
       fac(i)=fac(i-1)*float(i-1)
      enddo
      
      llx=mx
      lly=my
      do lx=1,nrx-1
       call slinf(xr(lx),xt,knotx,llx)
       call bvald(xr(lx),xt,mx,llx,mx,dx)
       do ly=1,nry-1
        call slinf(yr(ly),yt,knoty,lly)
        call bvald(yr(ly),yt,my,lly,my,dy)
        do ix=1,mx
         do iy=1,my
          deriv=0.
          do sx=1,mx
           do sy=1,my
            deriv=deriv+f(nx*(lly-my+sy-1)+llx-mx+sx)*
     &           dx(mx*(sx-1)+ix)*dy(my*(sy-1)+iy)
           enddo                !sy
          enddo                 !sx
          s(mx*(my*((nrx-1)*(ly-1)+lx-1)+my-iy)+mx-ix+1)=
     &         deriv/fac(ix)/fac(iy)
         enddo                  !iy
        enddo                   !ix
       enddo                    !ly
      enddo                     !lx
      
      return
      
      end
      
c************************************************************************
      
      subroutine pp2d(x,xt,xr,nx,nrx,mx,xx,knotx,lx,dfxydx,f,s,fxy,init,
     &     y,yt,yr,ny,nry,my,yy,knoty,ly,dfxydy)
c     Version: 30/03/95
      
c     interpolation 2D equivalente a sbsp2d mais on calcule avec les
c     polynomes par morceaux d'ou un gain de temps justifie a partir
c     de quelques nx*ny appels
      
c     interpolation dans le tableau f(nx,ny) par produit tensoriel
c     de splines polynomiales d'ordres mx{y}>1, au point
c     (xx,yy) : x(1).le.xx.le.x(nx), y(1).le.yy.le.y(ny), on aura aussi
c     xt(lx).le.xx.lt.xt(lx+1)       yt(ly).le.yy.lt.yt(lt+1)
c     en entree prendre lx{y}=mx{y}
c     au premier appel, init=.false. il y a initialisation :
c     calul de knotx=2mx+nx, formation du tableau xt(knotx). idem pour y
c     !! ATTENTION le tableau f(nx,ny) initial est remplace par le
c     tableau des coefficients des B-splines en xy!!
c     ce calcul etant fait selon de Boor p.342
      
c     ax et ay doivent etre dimensionnes au moins a   max(nx*nx,ny*ny)
c     l et m  "       "       "       "       "       max(nx,ny)
      
c     les tables de travail xr, yr sont conservees car on peut faire appel
c     sequentiellement a des tables diverses
      
c     Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c     Version: 04/01/95
      
      implicit none
      
      integer dimab,dimxy,dimaxy,dimqxy,dimdxy
      parameter (dimab=500, dimxy=6)
      parameter (dimaxy=dimab**2, dimqxy=dimxy+1, dimdxy=dimxy**2)
      
      integer nx,ny,knotx,knoty,lx,ly,mx,my,l(dimab),m(dimab),i,j,k,
     &     llx,lly,ix,iy,sx,sy,nrx,nry
      
      double precision x(*),y(*),f(*),xt(*),yt(*),s(*),xx,yy,fxy
      double precision dfxydx,dfxydy,
     &     ax(dimaxy),ay(dimaxy),qx(dimqxy),qy(dimqxy),dx(dimdxy),
     &     dy(dimdxy),pivot,d2fdxy,deriv,fac(dimxy),xx1,yy1,
     &     xr(*),yr(*)
      
      logical init

      if(init)goto 100
      
c     initialisations
      
      if(nx .gt. dimab .or. mx .gt. dimxy)goto 200
      if(ny .gt. dimab .or. my .gt. dimxy)goto 200
      
      do i=1,dimaxy
       ax(i)=0.
       ay(i)=0.
      end do
      
      lx=mx
      ly=my
      
c     formation des tableaux xt, yt et points de raccord
      
      call snoein(x,xt,nx,mx,knotx)
      call snoein(y,yt,ny,my,knoty)
      nrx=nx-mx+2
      nry=ny-my+2
      do ix=1,nrx
       xr(ix)=xt(ix+mx-1)
      enddo
      do iy=1,nry
       yr(iy)=yt(iy+my-1)
      enddo
      
c     ay <- matrice 1-D des interpolations en y
      
      do j=1,ny
       call slinf(y(j),yt,knoty,ly)
       call bval(y(j),yt,my,ly,qy)
       do k=1,my
        ay(ny*(ly-my+k-1)+j)=qy(k)
       enddo                    !k
      enddo                     !j
      
c     ax <- transposee de ay, inversion de ax et ay=f*ax
      
      call strans(ay,ny,ax)
      
      call minv(ax,ny,pivot,l,m)
      
      if(pivot .eq. 0.)goto 220
      call smatpr(f,nx,ny,ax,ny,ay)
      
c     ax <- matrice 1-D des interpolations en x
      
      do i=1,dimaxy
       ax(i)=0.
      end do
      do j=1,nx
       call slinf(x(j),xt,knotx,lx)
       call bval(x(j),xt,mx,lx,qx)
       do k=1,mx
        ax(nx*(lx-mx+k-1)+j)=qx(k)
       enddo                    !k
      enddo                     !j
      
c     ax=inverse de ax et f=ax*ay
      
      call minv(ax,nx,pivot,l,m)
      if(pivot .eq. 0.)goto 230
      call smatpr(ax,nx,nx,ay,ny,f)
      
c     calcul des coefficients des pp. : 12-21 de schumaker
      
c     de facon a optimiser l'algoritme de horner avec derivee
c     dans s, les coefficients du polynome (derives/fac(.)) sont classes
c     par puissance decroissante :
c     fy3x3, fy3x2, fy3x, fy3; fy2x3, fy2x2, fy2x, fy2;
c     fyx3, fyx2, fyx,fy; fx3, fx2, fx, f
      
c     s(ix,iy,lx,ly)=s(mx*(my*((nrx-1)*(ly-1)+lx-1)+iy-1)+ix)
c     il faut ajouter 1 a l'ordre de la derivation, ainsi
c     pour fy2x on aura ix=mx-2, iy=my-3 et, pour l'ordre de derivation i,j
c     correspondant aux indices k=i+1, l=j+1 pour bvald, on aura
c     ix=mx-k+1, iy=my-l+1
      
      fac(1)=1.                 !calcul de (i-1)!
      fac(2)=1.
      do i=3,dimxy
       fac(i)=fac(i-1)*float(i-1)
      enddo
      
      llx=mx
      lly=my
      do lx=1,nrx-1
       call slinf(xr(lx),xt,knotx,llx)
       call bvald(xr(lx),xt,mx,llx,mx,dx)
       do ly=1,nry-1
        call slinf(yr(ly),yt,knoty,lly)
        call bvald(yr(ly),yt,my,lly,my,dy)
        do ix=1,mx
         do iy=1,my
          deriv=0.
          do sx=1,mx
           do sy=1,my
            deriv=deriv+f(nx*(lly-my+sy-1)+llx-mx+sx)*
     &           dx(mx*(sx-1)+ix)*dy(my*(sy-1)+iy)
           enddo                !sy
          enddo                 !sx
          s(mx*(my*((nrx-1)*(ly-1)+lx-1)+my-iy)+mx-ix+1)=
     &         deriv/fac(ix)/fac(iy)
         enddo                  !iy
        enddo                   !ix
       enddo                    !ly
      enddo                     !lx
      
c     localisation de (xx,yy)
      
 100    if(xx .lt. x(1))then
       xx1=xr(1)
       lx=1
       write(*,*)'dans pp2d extrapolation pour x=',xx,'.lt. x(1)=',x(1)
      else if(xx .gt. x(nx))then
       xx1=xr(nrx)
       lx=nrx-1
       write(*,*)'dans pp2d extrapolation pour x=',xx,'.gt. x(nx)=',
     &      x(nx)
      else
       call slinf(xx,xr,nrx,lx)
       xx1=xx
      endif
      
      if(yy .lt. y(1))then
       yy1=yr(1)
       ly=1
       write(*,*)'dans pp2d extrapolation pour y=',yy,'.lt. y(1)=',y(1)
      elseif(yy .gt. y(ny))then
       yy1=yr(nry)
       ly=nry-1
       write(*,*)'dans pp2d extrapolation pour y=',yy,'.gt. y(ny)=',
     &      y(ny)
      else
       call slinf(yy,yr,nry,ly)
       yy1=yy
      endif
      
c     interpolation
      
      do iy=1,my                !lly premier indice pour y**(my-iy)
       lly=mx*(my*((nrx-1)*(ly-1)+lx-1)+iy-1)+1
       call horder(mx,lly,s,xx1-xr(lx),ax(iy),ay(iy))
      enddo                     !iy
      
      call horder(my,1,ax,yy1-yr(ly),fxy,dfxydy)
      call horder(my,1,ay,yy1-yr(ly),dfxydx,d2fdxy)
      
      return
      
 200  write(*,*)'dans pp2d regarder le commentaire sur les dimensions'
      write(*,201)nx,ny,mx,my
 201  format(1x,'nx=',i4,' ou ny=',i3,' ou mx=',i3,' ou my=', i3,
     &     ' est trop grand')
      stop
      
 220  write(*,*)'dans l inversion 1 de ax'
      stop
      
 230  write(*,*)'dans l inversion 2 de ax'
      stop
      
      end
      
c******************************************************************
      
      subroutine snoein(x,y,n,m,knot)
c     Version: 04/01/95
      
c     determine la sequence de noeuds de raccord y(knot)
c     pour une interpolation "optimale"
c     par B-splines d'ordre m sur la suite strictement
c     croissante x de n points de donnee, cf. de Boor p.219 formule (10)
c     aux limites le polynome d'interpolation s'appuie sur m points de donnee
      
c     Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
      
      implicit none
      
      integer m,knot,n,i,j
      
      double precision x(*),y(*),mm1,eps
c     double precision x(120),y(120),mm1,eps
      data eps/1.d-12/
      
c     verification de la stricte croissance de la suite des x
      
      do i=1,n-1
       if(x(i) .ge. x(i+1))then
        write(*,*)'dans snoein la suite des abscisses n est pas',
     &       ' strictement croissante en i=',i
        write(*,*)'nombre de points : ',n
        write(*,*)'abscisses x=',x(i)
        write(*,*)(x(knot),knot=1,n)
        pause
       endif
      enddo
      
c     pour l'interpolation spline il faut n.ge.m
      
      if(n .lt. m)then
       write(*,11)n,m
 11    format(1x,'dans snoein n=',i4,'.lt.',i3,'=m')
       pause 'stop...'
      endif
      
      mm1=m-1
      
c     formation du tableau des y
      
      knot=0
      
      do i=1,m
       knot=knot+1
       y(knot)=x(1)-eps
      enddo
      
      do i=1,n-m
       knot=knot+1
       y(knot)=0.
       do j=i+1,i+m-1
        y(knot)=y(knot)+x(j)
       enddo
       y(knot)=y(knot)/mm1
      enddo
      
      do i=1,m
       knot=knot+1
       y(knot)=x(n)+eps
      enddo
      
      return
      
      end
      
c********************************************************************
      
      subroutine slinf(x,y,n,l)
      
c     y(n) : suite croissante (non strictement) des points de table
c     recherche de l'indice l tel que :
c     y(l) .le. x .lt. y(l+1)
      
c     s'il y a debordement a gauche ou a droite on prend:
c     x < y(1) <= y(l) < y(l+1) ou y(l) < y(l+1) <= y(n) < x
      
c     on commence la recherche au voisinage de l
      
c     Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c     Version: 04/01/95
      
      implicit none
      
      integer n,l
      
      double precision y(*),x
      
      
      if(x .le. y(1))then
       l=1
       if(x .lt. y(1))write(*,*)'slinf x .lt. y(1)',x,' < ',y(1)
       do while(y(l+1) .le. y(1))
        l=l+1
       enddo
      elseif(x .ge. y(n))then
       l=n-1
       if(x .gt. y(n))write(*,*)'slinf x .gt. y(n)',x,' > ',y(n)
       do while(y(l) .ge. y(n))
        l=l-1
       enddo
      else
       l=max(1,min(l,n-1))      !initialisation de l
       if(x .ge. y(l))then
        do while(x .ge. y(l+1))
         l=l+1
        enddo
       else
        do while(x .lt. y(l))
         l=l-1
        enddo
       endif
      endif
      
      return
      
      end
      
      
c**********************************************************************
      
      subroutine bval(x,y,m,l,q)
      
c     calcul les B-splines normalisees d'ordre m>1
c     au point x tel que y(l).le.x.lt.y(l+1)
c     si x est l'extremite droite : y(l).lt.x.eq.y(l+1) on obtient
c     la limite a droite
c     les valeurs des m b-splines non nulles : l-m+1, ..., l sont 
c     dans q(1), ..., q(m)
c     !! ATTENTION : q doit etre dimensionne au moins a m+1 !!
      
c     d'apres l'algorithme 5-5 p.192 de schumaker
      
c     Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c     Version: 04/01/95
      
      implicit none
      
      integer*4 m,l
      integer*4 j,i
      
      double precision y(*),q(*),x
      double precision denom,a1,a2
      
      do j=1,m-1
       q(j)=0.
      end do
      
      q(m)=1./(y(l+1)-y(l))
      q(m+1)=0.
      
      do j=2,m-1
       do i=m-j+1,m
        denom=y(i+l-m+j)-y(i+l-m)
        a1=(x-y(i+l-m))/denom
        a2=1.d0-a1
        q(i)=a1*q(i)+a2*q(i+1)
       enddo
      enddo
      
      do i=1,m
       q(i)=(x-y(i+l-m))*q(i)+(y(i+l)-x)*q(i+1)
      enddo
      
      return
      
      end
      
c*********************************************************************
      
      subroutine strans(a,n,b)
      
c     b = transposee de a(n,n)
      
c     Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c     Version: 04/01/95
      
      implicit none
      
      integer*4 n
      integer*4 i,j
      
      double precision a(*),b(*)
      
      do i=1,n
       do j=1,n
        b(n*(j-1)+i)=a(n*(i-1)+j)
       enddo
      enddo
      
      return
      
      end
      
c**************************************************************
      
      subroutine minv(a,n,d,l,m)
      
c     inversion de la matrice a, programme de la ssp 360
      
c     Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c     Version: 04/01/95
c     Modifications: 25/05/91
      
      implicit none
      
      integer*4 l(*),m(*),n
      integer*4 k,nk,kk,j,iz,ij,ki,ji,i,ik,jq,jr,jp,jk,kj
      
      double precision a(*),d
      double precision biga,hold
c     
c     search for largest element
c     
      d=1.
      
      nk=-n
      do k=1,n
       nk=nk+n
       l(k)=k
       m(k)=k
       kk=nk+k
       biga=a(kk)
       do j=k,n
        iz=n*(j-1)
        do i=k,n
         ij=iz+i
         if((abs(biga)-abs(a(ij))).lt.0) then
          biga=a(ij)
          l(k)=i
          m(k)=j
         endif
        enddo
       enddo
c     
c     interchange rows
c     
       j=l(k)
       if((j-k).gt.0) then
        ki=k-n
        do i=1,n
         ki=ki+n
         hold=-a(ki)
         ji=ki-k+j
         a(ki)=a(ji)
         a(ji) =hold
        enddo
       endif
c     
c     interchange columns
c     
       i=m(k)
       if((i-k).gt.0) then
        jp=n*(i-1)
        do j=1,n
         jk=nk+j
         ji=jp+j
         hold=-a(jk)
         a(jk)=a(ji)
         a(ji) =hold
        enddo
       endif

c     divide column by minus pivot (value of pivot element is
c     contained in biga)
          
       if((biga).eq.0.) then
        d=0.
        write(*,1000)biga
 1000   format(1x,'matrice non inversible dans minv biga=',1pe10.3)
        return
       endif
          
       do i=1,n
        if((i-k).ne.0) then
         ik=nk+i
         a(ik)=a(ik)/(-biga)
        endif
       enddo
c     
c     reduce matrix
c     
       do i=1,n
        ik=nk+i
        hold=a(ik)
        ij=i-n
        do j=1,n
         ij=ij+n
         if((i-k).ne.0) then
          if((j-k).ne.0) then
           kj=ij-i+k
           a(ij)=hold*a(kj)+a(ij)
          endif
         endif
        enddo
       enddo
c     
c     divide row by pivot
c     
       kj=k-n
       do j=1,n
        kj=kj+n
        if((j-k).ne.0) a(kj)=a(kj)/biga
       enddo
c     
c     replace pivot by reciprocal
c     
       a(kk)=1./biga
      enddo
c     
c     final row and column interchange
c     
      k=n-1
      do while (k.gt.0)
       i=l(k)
       if((i-k).gt.0) then
        jq=n*(k-1)
        jr=n*(i-1)
        do j=1,n
         jk=jq+j
         hold=a(jk)
         ji=jr+j
         a(jk)=-a(ji)
         a(ji) =hold
        enddo
       endif
       j=m(k)
       if((j-k).gt.0) then
        ki=k-n
        do i=1,n
         ki=ki+n
         hold=a(ki)
         ji=ki-k+j
         a(ki)=-a(ji)
         a(ji)=hold
        enddo
       endif
       k=(k-1)
      enddo
            
      return      
      end
      
c*************************************************************************
      
      subroutine smatpr(a,la,calb,b,cb,ab)
      
c     produit de matrices ab=a*b
c     la : nombre de lignes de a
c     calb : nombre de colonnes de a et de lignes de b
c     cb : nombre de colonnes de b
      
c     Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c     Version: 04/01/95
      
      implicit none
      
      integer*4 la,calb,cb
      integer*4 i,j,k
      
      double precision a(*),b(*),ab(*)
      
      do i=1,la
       do j=1,cb
        ab(la*(j-1)+i)=0.
        do k=1,calb
         ab(la*(j-1)+i)=ab(la*(j-1)+i)+a(la*(k-1)+i)*b(calb*(j-1)+k)
        enddo
       enddo
      enddo
      
      return
      
      end
      
c********************************************************************
      
      subroutine bvald(x,y,m,l,r,d)
      
c     d(m*(j-1)+r):=d(r,j):= derivee d'ordre r-1,  1 .le. r .le. m, de
c     la l-m+j-ieme B-splines, 1 .le. j .le. m, d'ordre m non nulle
c     au point x, y(l) .le.  x .lt.  y(l+1).
c     les r derivees de 0 a r-1 sont calculees, et on a d(m,m)
      
c     Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c     Version: 04/01/95
      
      implicit none
      
      integer*4 m,l,r
      integer*4 j,i,mj,k,kml
      
      double precision x,y(*),d(*)
      double precision cd(10,0:10),q(11),n(10,10),a1,a2,denom
      
c     tests de dimension
      
      if(m .ge. 11)then
       write(*,*)'dans bvald, modifier les dimensions de cd(',
     &      m,'0:',m,'), de q(',m+1,') et de n(',m,',',m,')'
       stop
      endif
      
      if(r .gt. m)then
       write(*,*)'dans bvald, r=',r,'>m=',m
       stop
      endif
      
c     calcul des B-splines d'ordre 1 a m non nulles sur [y(l),y(l+1)[
c     adapte de schumaker 5-5
      
      do i=1,10
       do j=1,10
        n(i,j)=0.
       enddo
      enddo
      
      do i=1,m+1
       q(i)=0.
      enddo
      q(m)=1./(y(l+1)-y(l))
      do i=1,m
       n(1,i)=q(i)*(y(l-m+i+1)-y(l-m+i))
      enddo
      
      do j=2,m
       do i=m-j+1,m
        denom=y(i+l-m+j)-y(i+l-m)
        a1=(x-y(i+l-m))/denom
        a2=1.-a1
        q(i)=a1*q(i)+a2*q(i+1)
        n(j,i)=q(i)*(y(l-m+i+j)-y(l-m+i))
       enddo
      enddo
      
c     c(j,i) pour chacune des B-splines non nulle sur [y(l),y(l+1)[
c     adapte de schumaker 5-10
      
      do k=l-m+1,l
       kml=m*(k-l+m-1)
       do i=1,10
        do j=0,10
         cd(i,j)=0.
        enddo                   !j
       enddo                    !i
       cd(1,1)=1.
       
       do j=2,r
        mj=m-j+1
        do i=1,j
         denom=y(k+i-1+mj)-y(k+i-1)
         if(denom.eq.0.)then
          cd(j,i)=0.
         else
          cd(j,i)=float(mj)/denom*(cd(j-1,i)-cd(j-1,i-1))
         endif
        enddo                   !i
       enddo                    !j
       
c     do i=1,r
c     write(*,1000)(cd(i,j),j=1,r)
c 1000  format(1x,1p6d10.3)
c     enddo
c     write(*,*)' '
       
       do j=1,r
        d(kml+j)=0.
        do i=1,j
         d(kml+j)=d(kml+j)+cd(j,i)*n(m-j+1,k-l+m+i-1)
        enddo                   !i
       enddo                    !j
      enddo                     !k
      
      return
      
      end
      
c**********************************************************
      
      subroutine horder(m,l,a,x,p,d)
      
c     calcul de p(x) = a(l) * x**(m-1) + a(l+1) * x**(m-2) + ... + a(l+m-1)
c     et de d=dp/dx algorithme de Horner cf. burlirsh et stoer p. 270
      
c     Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c     Version: 04/01/95
      
      implicit none
      
      integer l,i,m
      
       double precision d,p,x,a(*)
      
      p=a(l)
      d=p
      do i=1,m-2
       p=p*x+a(l+i)
       d=d*x+p
      enddo                     !i
      p=p*x+a(l+m-1)
      
      return
      
      end



c**********************************************************************

      subroutine myclouds(p,t,sigma,a,kcloud)
	implicit double precision (a-h,o-z)
c     Version: 02/02/96
c Calculation of condensation in the atmosphere of Gl229
c These calculations are possible thanks to data provided by K. Lodders and B. Fegley
c Caution: All this is VERY rough!
c (c) Tristan Guillot
c Entrees:
c     p,t,sigma: pressure (dyn/cm2), temperature (K), wavenumber (cm-1)
c     a: assumed radius of cloud particles in microns (0.01<a<1)
c Sortie:
c     kcloud: scattering by clouds in cm2/g

       double precision kcloud

      parameter (nel=7)
	dimension fmol(nel),alph(nel),x(nel)
      character*20 cmol(nel)
      data cmol/'NH4H2PO4','ZnS','K2S','Na2S','MnS','Cr','MgSiO3'/
      data fmol/7d-7,9d-8,1d-7,2d-6,7d-7,1d-6,7d-5/
      data alph/115.0,126.0,127.9,128.0,134.4,135.8,136.6/

      xx=0.d0
      do i=1,nel
       x(i)=10**(-alph(i)+41*dlog10(t))/p*1d6
       if (x(i).gt.fmol(i)) x(i)=0.d0
       xx=xx+x(i)
      enddo
c      write(*,'(1p7d10.3)')(x(i),i=1,nel)

      sig2=sigma*sigma*1.d-16   !en Angstrom-2
c     Dalgarno & Williams, Proc.Phys.Soc. (1965)
      sray_h2=8.140d-13*sig2*sig2*(1+sig2*(1.572d6+sig2*(1.981d12
     &      +sig2*(2.307d18+sig2*(2.582d24+sig2*2.822d30)))))

      kcloud=6.022d23/2.2*xx*(a/2d-4)**3*sray_h2

      return
      end

c**********************************************************************

      subroutine opa_hmff(sigma,t,sff_hm)
c INPUT parameters:  
c sigma: wavenumber in cm-1
c t:  temperature in K
c
c OUTPUT parameter:
c  sff_hm: cross section in cm^5
c
c To get the opacity in cm^2/g, you then need to multiply sff_hm by N(H)*N(e)/rho
c where rho is the density in g/cm^3 and the N(i) are number densities in cm-3.  Note
c that the correction factor for stimulated emission is included in sff_hm.


c     Version: 24/4/02
c     Section efficace d'absorption free-free de H- d'apres Bell & Berrington (1987)
c     The factor for stimulated emission is included here
c     Auteur: D. Saumon
c Entree:
c     sigma: nombre d'onde [cm-1]
c     t: temperature en K
c Sortie:
c     sff_hm: coefficient d'absorption [cm2]

      implicit none
      double precision sigma,sff_hm,al,aj1(6),bj1(6),cj1(6),dj1(6)
      double precision ej1(6), fj1(6),aj2(6),bj2(6),cj2(6),dj2(6)
      double precision ej2(6),fj2(6),t,tcoeff, hj(6)
      double precision pi,echarg,kbol,hpl,clight,eve,amu,me,
     &       g,aradia,secon6,ah,ahe4,nl
      integer k,i

      DATA AJ1/ 0.D0, 2483.346, -3449.889, 2200.040, -696.271, 88.283/
      DATA BJ1/ 0.D0, 285.827, -1158.382, 2427.719, -1841.400, 444.517/
      DATA CJ1/ 0.D0, -2054.291, 8746.523, -13651.105, 8624.970,
     &          -1863.864/
      DATA DJ1/ 0.D0, 2827.776, -11485.632, 16755.524, -10051.530,
     &          2095.288/
      DATA EJ1/ 0.D0, -1341.537, 5303.609, -7510.494, 4400.067,
     &          -901.788/
      DATA FJ1/ 0.D0, 208.952, -812.939, 1132.738, -655.020, 132.985/
      DATA AJ2/ 518.1021, 473.2636, -482.2089, 115.5291, 2*0.D0/
      DATA BJ2/ -734.8666, 1443.4137, -737.1616, 169.6374, 2*0.D0/
      DATA CJ2/ 1021.1775, -1977.3395, 1096.8827, -245.649, 2*0.D0/
      DATA DJ2/ -479.0721, 922.3575, -521.1341, 114.243, 2*0.D0/
      DATA EJ2/ 93.1373, -178.9275, 101.7963, -21.9972, 2*0.D0/
      DATA FJ2/ -6.4285, 12.3600, -7.0571, 1.5097, 2*0.D0/

C  AL is the wavelength in microns.  The fit can be extended to ~ 20 microns and down to T ~ 800K
C  without serious difficulty.
      AL=1.d4/sigma
      sff_hm=0.D0
      if(t.lt.800 .or. al.gt.20) return
      TCOEFF=5040.0/T
      DO K=1,6
        IF(AL.GT.0.3645) THEN
            HJ(K)=1.D-29*(AL*AL*AJ1(K) + BJ1(K) + (CJ1(K) + (DJ1(K) +
     !               (EJ1(K) + FJ1(K)/AL)/AL)/AL)/AL)
          ELSE
            IF(AL.LT.0.1823) AL=0.1823D0
            HJ(K)=1.D-29*(AL*AL*AJ2(K) + BJ2(K) + (CJ2(K) + (DJ2(K) +
     !               (EJ2(K) + FJ2(K)/AL)/AL)/AL)/AL)
        ENDIF
      ENDDO   
      DO I=1,6
        sff_hm=sff_hm + TCOEFF**((I+1.)/2.)*HJ(I)
      ENDDO
	  kbol = 1.380658d-16
      sff_hm=sff_hm*kbol*t
      return
      end

