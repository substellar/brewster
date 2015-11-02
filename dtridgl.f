      SUBROUTINE DTRIDGL(L,AF,BF,CF,DF,XK)
C DOUBLE PRESCISION VERSION OF TRIDGL
      IMPLICIT none 
      integer nmax, i,l
      PARAMETER (NMAX=301)
      double precision AF(L),BF(L),CF(L),DF(L),XK(L)
      double precision AS(NMAX),DS(NMAX)
      double precision x,xkb
C* THIS SUBROUTINE SOLVES A SYSTEM OF TRIDIAGIONAL MATRIX
C*  EQUATIONS. THE FORM OF THE EQUATIONS ARE:
C*  A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = D(I)
C* WHERE I=1,L  LESS THAN 103.
C* ..............REVIEWED -CP........
      AS(L) = AF(L)/BF(L)
      DS(L) = DF(L)/BF(L)
      DO 10 I=2,L
           X=1./(BF(L+1-I) - CF(L+1-I)*AS(L+2-I))
           AS(L+1-I)=AF(L+1-I)*X
           DS(L+1-I)=(DF(L+1-I)-CF(L+1-I)*DS(L+2-I))*X
   10 CONTINUE
      XK(1)=DS(1)
      DO 20 I=2,L
           XKB=XK(I-1)
           XK(I)=DS(I)-AS(I)*XKB
   20 CONTINUE
      RETURN
      END

