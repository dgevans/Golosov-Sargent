C***************************************************************************
C   See the MATLAB file doubles.m for comments on use
C
C   CODE LAST MODIFIED:  Dec 27, 1998
C   DOCUMENTATION LAST MODIFIED:  September 7, 1996
C
C***************************************************************************


#include "lib1.f"
#include "lib2.f"

#define G_IN   PRHS(1)
#define D_IN   PRHS(2)
#define H_IN   PRHS(3)
#define TOL_IN   PRHS(4)
#define MAXIT_IN PRHS(5)

#define V_OUT PLHS(1)
#define INFO_OUT PLHS(2)

#define NUMIN 5
#define NUMOUT 2


C**********************************************************************
C
C     GATEWAY FUNCTION FOR THE DOUBLING ALGORITHM FOR SYLVESTER
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
      INTEGER*4 M,N
      INTEGER*4 WORK1,WORK2,WORK3,WORK4
c-----------------------------------------------------------------------
      INTEGER*4 MAXIT 
	

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
      ELSEIF ((1 .NE. MXGETM(MAXIT_IN)) .OR.
     x         (1 .NE. MXGETN(MAXIT_IN))) THEN
         CALL MEXERRMSGTXT ('Fifth input must be a scalar')
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
C  Create workspace varaibles
      WORK1 = MXCALLOC(M*M,8)
      WORK2 = MXCALLOC(N*N,8)
      WORK3 = MXCALLOC(M*N,8)
      WORK4 = MXCALLOC(M*N,8)
C
C put the inputs in the proper from 
      G = MXCALLOC(M*M,8)
      H = MXCALLOC(N*N,8)
      CALL MCOPYM(%VAL(MXGETPR(G_IN)),%VAL(G),M,M)
      CALL MCOPYM(%VAL(MXGETPR(H_IN)),%VAL(H),N,N)
      TOL = MXGETPR(TOL_IN)
      MAXIT = MXGETPR(MAXIT_IN)
c
c call the computational routine
      CALL DOUBCJ2(%VAL(G),%VAL(V),%VAL(H),
     x          %VAL(TOL),%VAL(MAXIT),M,N,
     x          %VAL(WORK1),%VAL(WORK2),%VAL(WORK3),%VAL(WORK4))
C

C
C Finish
      RETURN
      END

C*****************************************************************
C
C   computational routine for the doubling algorithm for Sylvester equations
C
C    PARAMETERS:
C         A      MxM double precision matrix
C         V      (On Input) MxN double precision matrix
C                (On Output) contains the solution
C         D      NxN double precision matrix
C         TOL    double precision:  criteria for convergence
C         MAXIT  double precision: maximum number of iterations
C         M      integer*4
C         N      integer*4   
C         WORKI  N dimesnional integer*4 vector
C         WORK   NxN double precision matrix
C         WORK2  NxN double precision matrix
C         WORK3  NxN double precision matrix
C         WORK4  NxN double precision matrix
C       
C*****************************************************************



      SUBROUTINE  DOUBCJ2(A,V,D,TOL,MAXIT,M,N,WORK1,WORK2,WORK3,WORK4)

      DOUBLE PRECISION   A(M,M), V(M,N), D(N,N),TOL
      DOUBLE PRECISION   WORK1(M,M),WORK2(N,N),WORK3(M,N),WORK4(M,N)
      DOUBLE PRECISION   MAXIT
      INTEGER*4  M,N

      LOGICAL MCRITR
      INTEGER*4  ITER

      ITER=1
C Compute WORK4 = A*V*D
   10 CALL MMUL (V,D,WORK3,M,N,N)
      CALL MMUL (A,WORK3,WORK4,M,M,N)

C Check to see if finished
      CALL MADD(V,WORK4,M,N)
      IF (MCRITR(WORK4,V,TOL,M,N)) GO TO 100
      IF (ITER .GE. MAXIT) GO TO 100
      ITER = ITER + 1

C  Compute a*a  and d*d

      CALL MMUL(A,A,WORK1,M,M,M)
      CALL MMUL(D,D,WORK2,N,N,N)

C Compute WORK4 = A*V*D
      CALL MMUL (V,WORK2,WORK3,M,N,N)
      CALL MMUL (WORK1,WORK3,WORK4,M,M,N)


C Check to see if finished
      CALL MADD(V,WORK4,M,N)
      IF (MCRITR(WORK4,V,TOL,M,N)) GO TO 100
      IF (ITER .GE. MAXIT) GO TO 100
      ITER = ITER + 1 

C  Compute a*a  and d*d

      CALL MMUL(WORK1,WORK1,A,M,M,M)
      CALL MMUL(WORK2,WORK2,D,N,N,N)

C Loop back
      GO TO 10


C Finish
100   CONTINUE
      END

