C***************************************************************************
C   See the MATLAB file sylvest.m for comments on use
C
C   CODE LAST MODIFIED:  December 27, 1998
C   DOCUMENTATION LAST MODIFIED:  September 7, 1996
C
C***************************************************************************

 
#include "dblas2.f"
#include "lib1.f"
#include "lib2.f"
#include "lib3.f"
#include "lib4.f"

#define G_IN   PRHS(1)
#define D_IN   PRHS(2)
#define H_IN   PRHS(3)
#define TOL_IN   PRHS(4)

#define V_OUT PLHS(1)
#define INFO_OUT PLHS(2)

#define NUMIN 4
#define NUMOUT 2


C**********************************************************************
C
C     GATEWAY FUNCTION FOR THE HESSENBERG-SCHUR FOR SYLVESTER
C     EQUATIONS
C
C*********************************************************************


      SUBROUTINE MEXFUNCTION(NLHS, PLHS, NRHS, PRHS)
      INTEGER*4 PLHS(*), PRHS(*)
      INTEGER*4 NLHS, NRHS
      INTEGER*4 MXCREATEFULL, MXGETPR
      INTEGER*4 MXGETM, MXGETN
C-----------------------------------------------------------------------
      INTEGER*4 G,H,TOL
      INTEGER*4 V,INFO
C----------------------------------------------------------------------
      INTEGER* 4 M,N,S,M2
      INTEGER*4 GHESS,A,A2,TEMP,WORKI

C Check the number of inputs and outputs
      IF (NRHS .NE. NUMIN) THEN
        CALL MEXERRMSGTXT('Wrong number of input parameters')
      ELSEIF (NLHS .NE. NUMOUT) THEN
        CALL MEXERRMSGTXT('Wrong number of output parameters')
      ENDIF
C
C Check the sizes of the inputs
      M = MXGETM(D_IN)
          N = MXGETN(D_IN)
      IF  ( (M  .NE.  MXGETM(G_IN))  .OR.
     x      (M  .NE.  MXGETN(G_IN))) THEN
        CALL MEXERRMSGTXT ('Bad dimensions for first input')
      ELSEIF ( (N .NE.  MXGETM(H_IN))  .OR.
     x         (N .NE.  MXGETN(H_IN))) THEN
        CALL MEXERRMSGTXT ('Bad dimensions for third input')
      ELSEIF  ((1 .NE. MXGETM(TOL_IN)) .OR.
     x          (1 .NE. MXGETN(TOL_IN))) THEN
         CALL MEXERRMSGTXT ('Fourth input must be a scalar')
      ELSEIF ( M .LT. N) THEN
             CALL MEXERRMSGTXT ('Inputs of the wrong order')
      ENDIF
C
C Check the inputs for some common errors
      CALL CHECKIN (NUMIN,PRHS)
C
C  Initialize the outputs
      V_OUT = MXCREATEFULL(M,N,0)
      INFO_OUT = MXCREATEFULL(1,1,0)
C
      V = MXGETPR(V_OUT)
      INFO = MXGETPR(INFO_OUT)
      CALL MCOPYM (%VAL(MXGETPR(D_IN)),%VAL(V),M,N)
C
C  Create variables for the workspace S could be S = S-1, that
C  is dim(A) = M(M+1)/2 + M - 1 
      M2 = 2*M
      S = (M*(M+1)/2) + M
      GHESS = MXCALLOC(S,8)
      A = MXCALLOC(S,8)
      TEMP = MXCALLOC(M*N,8)
      A2 = MXCALLOC(M2*M2,8)
      WORKI = MXCALLOC(M2,4)
C
C put the inputs in the proper from . Note that G,H and TOL are not 
C overwritten in the computational section
      G = MXGETPR(G_IN)
      H = MXGETPR(H_IN)
      TOL = MXGETPR(TOL_IN)
      CALL FULL2HESS(%VAL(G),%VAL(GHESS),M)
c
c call the computational routine
      CALL SYLVEST (%VAL(G),%VAL(V),%VAL(H),M,N,%VAL(GHESS),
     x              %VAL(A),%VAL(TEMP),%VAL(A2),%VAL(WORKI),
     x              %VAL(TOL),%VAL(INFO))
C

C finish
      RETURN
      END


C*******************************************************************
C
C     Computational routine for solving a sylvester equation with
C     a Hessenberg-Schur algorithm.ine for solving a sylvester equation
C
C********************************************************************

     

      SUBROUTINE SYLVEST(G,D,H,M,N,GHESS,A,TEMP,A2,
     x                   WORKI, TOL, INFO)
      DOUBLE PRECISION G(M,M), H(N,N), D(M,N),TOL,INFO
      DOUBLE PRECISION GHESS(*),A(*), TEMP(M,N), A2(2*M,2*M)
      INTEGER*4  M,N,WORKI(*)
          
      INTEGER*4 I,L,U
      LOGICAL  FAIL
          
      INFO = 0.0D0
c     I = 1
c
c  handle I = 1 outside
      IF (N .EQ. 1)  GOTO 95
      I = 1
      L = 0
      IF (DABS(H(I+1,I)) .LT. TOL) THEN
         CALL FORMA(A,GHESS,H(I,I),M)
         CALL VLINEARH(A,D(1,I),M,WORKI,FAIL)
         IF (FAIL) GOTO 100
         I = I + 1
      ELSE
        CALL FORMA2(A2,G,H(I,I),H(I+1,I),H(I,I+1),H(I+1,I+1),M)
        CALL VLINEAR(A2,D(1,I),2*M,WORKI,FAIL)
        IF (FAIL) GOTO 100
        I = I + 2
      ENDIF
c
c  handle the rest of the I's with the possible exception of I = N
      IF (I .GE. N) GO TO 95
   
c
c  compute temp  check setting of L
   10 U = I-1
      CALL MMUL(G,D(1,L+1),TEMP(1,L+1),M,M,U-L)
      L = U
      IF (DABS(H(I+1,I)) .LT. TOL) THEN
          CALL FORMA(A,GHESS,H(I,I),M)
          CALL MMULADDW (TEMP,H(1,I),D(1,I),M,U,1,M,N,M)
          CALL VLINEARH(A,D(1,I),M,WORKI,FAIL)
          IF (FAIL) GOTO 100
          I = I + 1       
      ELSE
          CALL FORMA2(A2,G,H(I,I),H(I+1,I),H(I,I+1),H(I+1,I+1),M)
          CALL MMULADDW (TEMP,H(1,I),D(1,I),M,U,1,M,N,M)
          CALL MMULADDW (TEMP,H(1,I+1),D(1,I+1),M,U,1,M,N,M)
          CALL VLINEAR(A2,D(1,I),2*M,WORKI,FAIL)
          IF (FAIL) GOTO 100
          I = I + 2
      ENDIF      
   90 IF (I .LT. N) GO TO 10
           
c  
c handle I = N if applicable
   95 IF (I .EQ. N) THEN
        U = I-1
        CALL MMUL(G,D(1,L+1),TEMP(1,L+1),M,M,U-L)
        CALL FORMA(A,GHESS,H(I,I),M)
        CALL MMULADDW (TEMP,H(1,I),D(1,I),M,U,1,M,N,M)
        CALL VLINEARH(A,D(1,I),M,WORKI,FAIL)
        IF (FAIL) GOTO 100               
      ENDIF
      RETURN
c
c handle a singular A matrix
  100 INFO = REAL(I)

      RETURN
      END

                 
C*************************************************************************
C
C     SUBROUTINE:  FORMA
C
C     SYNOPSIS:  Forms the A matrix in SYLVEST
C
C******************************************************************************
                 
      SUBROUTINE FORMA(A,G,H1,N)
      DOUBLE PRECISION A(*),G(*),H1
      INTEGER*4 N
          
      INTEGER*4  I,K,S
c
c  multiply the matrix G by (-1)*H1
      S = (N*(N+1)/2) + N - 1
      DO 10 I =1,S
        A(I) = - G(I)*H1
   10 CONTINUE
c
c  Add the identity matrix to the result
      K = -1
      DO 20 I = 1,N
         K = K+I+1   
         A(K) = A(K) + 1
   20 CONTINUE
c
c  finish
      RETURN
      END
              



C*******************************************************************
C
C     SUBROUTINE: FORMA2
C
C     SYNOPSIS:  Forms the A2 matrix for solving a sylvester equation
C
C**********************************************************************

      SUBROUTINE FORMA2(A2,G,H1,H2,H3,H4,M)
      DOUBLE PRECISION A2(2*M,2*M),G(M,M),H1,H2,H3,H4
      INTEGER*4  M
          
      DOUBLE PRECISION X
      INTEGER*4  I,J
          
      DO 20 I=1,M
        DO 10 J=1,M
        X = - G(I,J)
        A2(I,J) = X*H1
        A2(I,J+M) = X*H2
        A2(I+M,J) = X*H3
        A2(I+M,J+M) = X*H4
   10   CONTINUE
   20 CONTINUE
   
      CALL MINCI(A2,2*M)
          
      RETURN
      END
          
