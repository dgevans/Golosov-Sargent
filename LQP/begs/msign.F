C***************************************************************************
C   See the MATLAB file msign.m for comments on use
C
C   CODE LAST MODIFIED:  December 27, 1998
C   DOCUMENTATION LAST MODIFIED:  September 7, 1996
C
C***************************************************************************




#include "dblas2.f"
#include "lib1.f"
#include "lib2.f"
#include "lib3.f"
 


#define H_IN  PRHS(1)
#define TOL_IN PRHS(2)
#define MAXIT_IN PRHS(3)

#define Z_OUT   PLHS(1)
#define INFO_OUT  PLHS(2)

#define NUMIN 3
#define NUMOUT 2

C***********************************************************************
C
C   Gateway function for calculation of the matrix sign
C
C***********************************************************************


      SUBROUTINE MEXFUNCTION(NLHS, PLHS, NRHS, PRHS)
      INTEGER*4 PLHS(*), PRHS(*)
      INTEGER*4 NLHS, NRHS
      INTEGER*4 MXCREATEFULL, MXGETPR
      INTEGER*4 MXGETM, MXGETN

      INTEGER*4  H,TOL,MAXIT
      INTEGER*4  WORK,WORK2,WORKI,WORKR
      INTEGER*4  Z, INFO

      INTEGER*4 N 

C
C  Check the number of inputs and outputs
      IF (NRHS .NE. NUMIN) THEN
        CALL MEXERRMSGTXT('Wrong number of input arguments')
      ELSEIF (NLHS .NE. NUMOUT) THEN
        CALL MEXERRMSGTXT('Wrong number of output arguments')
      ENDIF
C
C  Find out sizes of some of the inputs
      N = MXGETM(PRHS(1))
C
C  Check the sizes of the inputs
      IF (N .NE. MXGETN(H_IN)) THEN
       CALL MEXERRMSGTXT('First input must be square')
      ELSEIF ((MXGETM(TOL_IN) .NE. 1) .OR.
     x        (MXGETN(TOL_IN) .NE. 1)) THEN
       CALL MEXERRMSGTXT ('Second input must be a scalar')
      ELSEIF ((MXGETM(MAXIT_IN) .NE. 1) .OR.
     x        (MXGETN(MAXIT_IN) .NE. 1)) THEN
       CALL MEXERRMSGTXT ('Third input must be a scalar')
      ENDIF
C
C  Check the inputs for some common errors
      CALL CHECKIN (NUMIN,PRHS)
C
C Handle the inputs
      H = MXGETPR(H_IN)
      TOL = MXGETPR(TOL_IN)
      MAXIT = MXGETPR(MAXIT_IN)
C
C Handle the outputs
      Z_OUT = MXCREATEFULL(N,N,0)
      INFO_OUT = MXCREATEFULL(2,1,0)
      Z = MXGETPR(PLHS(1))
      INFO = MXGETPR(PLHS(2))
C
C  Allocate temporary work space.  
      WORKI = MXCALLOC (N,4)
      WORKR = MXCALLOC (N,8)
      WORK = MXCALLOC(N*N,8)
      WORK2 = MXCALLOC(N*N,8)
C
C  Copy inputs that doubc will modify
      CALL MCOPYS(%VAL(H),%VAL(Z),N)

C Call the computational routine
      CALL MSIGN (%VAL(Z), %VAL(TOL), %VAL(MAXIT),
     x             N,%VAL(WORK),%VAL(WORK2),%VAL(WORKR),
     x            %VAL(WORKI),%VAL(INFO))
C
C Bring the symmetric storage matrix Z up to full storage
      CALL MCONSR(%VAL(Z),N)
C
C Finish
      RETURN
      END

C_____________________________________________________________
C
c   Computational routine for the  matrix sign algorithm
C
C_________________________________________________________________


      SUBROUTINE MSIGN (Z,TOL,MAXIT,N,ZINV,WORK2,WORKR,
     x                  WORKI,INFO)
      DOUBLE PRECISION Z(N,N)
      DOUBLE PRECISION ZINV(N,N),WORK2(N,N), WORKR(N)
      DOUBLE PRECISION TOL, MAXIT
      INTEGER*4 N 
      DOUBLE PRECISION INFO(2)
      INTEGER*4 WORKI(N)
C
C     local variables
C    
      DOUBLE PRECISION  DETZ, INVN, C
      DOUBLE PRECISION  I
      LOGICAL ERR, DONE



C  Initialize variables for the main loop (includes copying the matrix Z)
      INFO(1) = 0
      INFO(2) = 0
      INVN = 1.0D0/N   
      I = 1    
      CALL MCOPYS(Z,ZINV,N)

C  Begin the main loop  -- Calculate the determinant and the inv of Z
   10 CALL MDETINVS(ZINV,N,DETZ,WORKR,WORKI,ERR)
      IF (ERR) GO TO 100
C     
C  Calculate c 
      C =DABS(DETZ) ** INVN  
C
C  Compute WORK2 = J*ZINV*J  
      CALL JTRANS(ZINV,WORK2,N)
C
C Find out if the maximum number of allowable iterations has been 
C performed
      IF (I .GE. MAXIT) THEN
         GO TO 120
      ENDIF 
C
C Find out if the new value of Z is within tolerance of the old value
C of Z and add the new part of Z to the old Z
      CALL DONESIGN(Z,WORK2,ZINV,(1/(2*C)),(C/2),N,TOL,DONE)
       I = I + 1
      IF (.NOT. DONE) GO TO 10
      GO TO 120 
C           
C Handle a error in matrix inversion
  100 INFO(2) = 1D0  
  120 INFO(1) = I 
      RETURN
      END


      SUBROUTINE JTRANS(A,B,N)
      INTEGER*4  N
      DOUBLE PRECISION  A(N,N), B(N,N)
C
C     Subroutine performs the matrix multiplication
C                 B = J*A*J  
C     where   J=[0 I;-I 0]
C
      INTEGER*4  I,J,S
C
C
      S = N/2
       DO 20 J=1,S
        DO 10 I=1,J
          B(I,J) = -A(I+S,J+S)
          B(I+S,J+S) = -A(I,J)
  10   CONTINUE 
  20  CONTINUE           

      DO 40 I = 1,S
       DO 30 J=S+1,N
         B(I,J) = A(J-S,I+S)
  30   CONTINUE
  40  CONTINUE
         
      RETURN
      END 






 
      SUBROUTINE DONESIGN(A,B,C,ASCALE,BSCALE,N,TOL,DONE)
      INTEGER*4  N
      DOUBLE PRECISION  A(N,N),B(N,N),C(N,N),ASCALE,BSCALE,TOL
      LOGICAL DONE
C
C    This subroutine performs the following computations
C           B = (ASCALE-1) *A + BSCALE*B
C           A = ASCALE*A + BSCALE*B
C           C = A
C    when it is assumed that the inputs A and B are both square and symmetric.
C    (This assumption is not checked.)
C
C     Upon exit DONE=.true. if  norm(B) is less than the norm(A)*tol. Otherwise
C     DONE=.false.
C
C   Parameters
C
C      A       NxN DOUBLE PRECISION MATRIX
C      B       NxN DOUBLE PRECISION MATRIX
C      ASCALE  Scalar double precision
C      BSCALE  Scalar double precision
C      N       Order of A and B -- integer
C      TOL     convergence criterion
C      DONE    boolean
C
C   Local Declarations
      INTEGER*4 I,J
      DOUBLE PRECISION X,RTOL
      DOUBLE PRECISION MNORM2
      LOGICAL MCRIT2
C
C  calculate (ASCALE-1)*A + BSCALE*B 
      X = ASCALE - 1
       DO 20 J=1,N
        DO 10 I =1,J
          B(I,J) = X*A(I,J) + BSCALE*B(I,J)
          A(I,J) = B(I,J) + A(I,J)
          C(I,J) = A(I,J)
   10   CONTINUE   
   20 CONTINUE         
c
c determine the relative tolerance and check convergence
      RTOL = MNORM2(A,N,N) * TOL 
      DONE = MCRIT2(B,RTOL,N,N)
C
C finish
      RETURN
      END

