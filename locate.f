C
C  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
	SUBROUTINE LOCATE(XX,N,X,J)
        implicit double precision (a-h,o-z)
C	Table searching routine from Numerical Recipes.  For
C       N=14 it is about twice as fast as the previous method.
	DIMENSION XX(N)
	
	JL=0
	JU=N+1
10	IF (JU-JL.GT.1) THEN
   	  JM=(JU+JL)/2
	  IF ((XX(N).GT.XX(1)).EQV.(X.GT.XX(JM))) THEN
	    JL=JM
	  ELSE
	    JU=JM
	  ENDIF
        GOTO 10
	ENDIF
	J=JL
	RETURN
	END
