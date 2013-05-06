C***************************************************************************
C   See the MATLAB file doubler.m for comments on use
C
C   CODE LAST MODIFIED:  June 1, 1995
C   DOCUMENTATION LAST MODIFIED:  September 7, 1996
C
C***************************************************************************



#include "dblas2.f"
#include "lib1.f"
#include "lib2.f"
#include "lib3.f"


#define A_IN PRHS(1)
#define B_IN PRHS(2)
#define V_IN PRHS(3)
#define TOL_IN PRHS(4)
#define MAXIT_IN PRHS(5)

#define V_OUT PLHS(1)
#define INFO_OUT  PLHS(2)

#define NUMIN  5
#define NUMOUT  2


C**********************************************************************
C
C     GATEWAY FUNCTION FOR THE DOUBLING ALGORITHM   
C
C*********************************************************************

      SUBROUTINE MEXFUNCTION(NLHS, PLHS, NRHS, PRHS)
      INTEGER*4 PLHS(*), PRHS(*)
      INTEGER*4 NLHS, NRHS
      INTEGER*4 MXCREATEFULL, MXGETPR
      INTEGER*4 MXGETM, MXGETN

      INTEGER*4 N,S
      INTEGER*4 V,INFO,A,B,TOL
      INTEGER*4 WORKI,WORK,WORK2,WORK3,WORK4,WORK5 

C
C Check the number of inputs and outputs
      IF (NRHS .NE. NUMIN) THEN
        CALL MEXERRMSGTXT('Wrong number of input parameters')
      ELSEIF (NLHS .NE. NUMOUT) THEN
        CALL MEXERRMSGTXT('Wrong number of output parameters')
      ENDIF
C
C Check the sizes of the inputs
      N = MXGETM(A_IN)
      IF (N .NE. MXGETN(A_IN)) THEN
       CALL MEXERRMSGTXT('First input must be square')
      ELSEIF ((N .NE. MXGETM(B_IN))  .OR.
     x        (N .NE. MXGETN(B_IN)))  THEN
       CALL MEXERRMSGTXT('1st and 2nd inputs must be of the same order')
      ELSEIF ((N .NE. MXGETM(V_IN))  .OR.
     x        (N .NE. MXGETN(V_IN)))  THEN
       CALL MEXERRMSGTXT('1st and 3rd inputs must be of the same order')
      ELSEIF ((MXGETM(TOL_IN) .NE. 1) .OR.
     x        (MXGETN(TOL_IN) .NE. 1)) THEN
       CALL MEXERRMSGTXT ('Fourth input must be a scalar')
      ELSEIF ((MXGETM(MAXIT_IN) .NE. 1) .OR.
     x        (MXGETN(MAXIT_IN) .NE. 1)) THEN
       CALL MEXERRMSGTXT ('Fifth input must be a scalar')
      ENDIF
C
C Check the inputs for some common errors
      CALL CHECKIN (NUMIN,PRHS)
C
C Handle the Outputs            
      V_OUT = MXCREATEFULL(N,N,0)
      INFO_OUT = MXCREATEFULL(2,1,0)
      V = MXGETPR(V_OUT)
      INFO = MXGETPR(INFO_OUT)      
C
C Handle the inputs
      S = N*N
      A = MXCALLOC(S,8)
      B = MXCALLOC(S,8) 
      CALL MCOPYM(%VAL(MXGETPR(A_IN)),%VAL(A),N,N)
      CALL MCOPYM(%VAL(MXGETPR(B_IN)),%VAL(B),N,N)
      CALL MCOPYM(%VAL(MXGETPR(V_IN)),%VAL(V),N,N)
      TOL = MXGETPR(TOL_IN)
C
C  Allocate temporary work space. 
      WORKI = MXCALLOC(N,4)
      WORK = MXCALLOC(S,8)
      WORK2 = MXCALLOC(S,8)
      WORK3 = MXCALLOC(S,8)
      WORK4 = MXCALLOC(S,8)
      WORK5 = MXCALLOC(S,8)
C
C Call the computational routine
      CALL DOUBC(%VAL(A),%VAL(B),%VAL(V),%VAL(TOL),%VAL(INFO),N,
     x           %VAL(WORKI),%VAL(WORK),%VAL(WORK2),%VAL(WORK3),
     x           %VAL(WORK4),%VAL(WORK5),
     x           %VAL(MXGETPR( MAXIT_IN))) 

C
C Finish
      RETURN
      END

C*****************************************************************
C
C    COMPUTATIONAL ROUTINE FOR THE DOUBLING ALGORITHM
C       information inside parenthesis
C       relate the variables to notation
C       in Anderson, Hansen, McGrattan and Sargent
C       "Mechanics of Forming and Estimating Dynamic
C       Linear Economies"
C
C    PARAMETERS:
C         A      NxN double precision matrix(On input
C                   this must be $\alpha_0$)
C         B      NxN double precision matrix(On input 
C                   this must be $\beta_0$)
C         V      ON INPUT: NxN double precision matrix ($\gamma_0$)
C                ON OUTPUT:  matrix in the value function(The
C                       limit of the doubling algorithm.)
C         TOL    double precision:  criteria for convergence
C                       ($\epsilon$)
C         INFO   2 dimensional double precision vector which on
C                    output gives information on failure.
C                    See doubler.m for a more complete 
C                    description of INFO
C         N      integer*4   dimension of the problem
C         WORKI  N dimesnional integer*4 vector
C         WORK   NxN double precision matrix
C         WORK2  NxN double precision matrix
C         WORK3  NxN double precision matrix
C         WORK4  NxN double precision matrix
C         WORK5  NxN double precision matrix 
C         MAXIT  double precision: maximum number of iterations
C
C*****************************************************************


      SUBROUTINE DOUBC(A,B,V,TOL,INFO,N,WORKI,WORK,WORK2,WORK3,
     x                 WORK4,WORK5,MAXIT)
      DOUBLE PRECISION   A(N,N), B(N,N), V(N,N),TOL,INFO(2)
      DOUBLE PRECISION   WORK(N,N),WORK2(N,N),WORK3(N,N)
      DOUBLE PRECISION   WORK4(N,N),WORK5(N,N),MAXIT
      INTEGER*4 N,WORKI(N)
C
C  Local Variables
      INTEGER*4  ITER
      LOGICAL DONE, FAIL
C
C  Initialize variables
      ITER = 1
      INFO(1) = 0       
      INFO(2) = 0
C
C Compute  IDENT+b*v    t = inv(Ident+b*v).  Here t=WORK
   10 CALL MMUL(B,V,WORK,N,N,N)
      CALL MINCI(WORK,N)

C
C Compute  WORK2 = t*a   
C Now compute WORK2 = t*a  WORK3 = t*b

       CALL MCOPYM(A,WORK2,N,N)
       CALL MCOPYM(B,WORK3,N,N)
       CALL MLINEAR2(WORK,WORK2,WORK3,N,N,WORKI,FAIL)
       IF (FAIL) THEN
          INFO(2) = 1
          RETURN
       ENDIF
C
C Compute v= v+a'*v*t*a  and see if the relative norm difference between
C vlhs and vrhs is greater than TOL .  WORK3, WORK4 are used in this
C section but are free aftewards

      CALL MMUL2 (A,V,WORK,N,N,N)
      CALL MMULMMSF(WORK,WORK2,WORK4,N,N) 
      CALL DONEDOUB2(V,WORK4,N,TOL,DONE)
      IF (DONE)  GOTO 100
      IF (ITER.GE. MAXIT) THEN
         GOTO 100
      ENDIF
      ITER = ITER +1
C
C Compute b=a*t*b*a'. WORK3, WORK4,WORK5 are used in this section but are free
C afterwards.

      CALL MMUL3(WORK3,A,WORK4,N,N,N)
      CALL MMULMMSF(A,WORK4,WORK5,N,N)
      CALL MADD (B,WORK5,N,N)
C
C Compute a = a*t*a

      CALL MMUL(A,WORK2,WORK3,N,N,N)
      CALL MCOPYM(WORK3,A,N,N)
      GO TO 10
C Finish
  100 INFO(1) = ITER 
      CONTINUE
      END

C*********************************************************************
C
C    SUBROUTINE:   DONEDOUB2
C
C    SYNOPSIS:  Computes the matrix addition A = A + B when it is assumed
C                that both A and B are square and symmetric. In addition
C                if norm(B) is less than  norm(A)*tol then
C                DONE= .TRUE. on exit.  Otherwise DONE=.FALSE.
C                The assumption that A and B are
C                symmetric is not checked.
C
C    PARAMETERS
C        A   NxN  double precision matrix
C        B   NxN  double precision matrix
C        N   integer*4   the dimension of the problem
C        TOL double precision
C        DONE  Boolean
C
C**************************************************************************

      SUBROUTINE DONEDOUB2(A,B,N,TOL,DONE)
      INTEGER*4  N
      DOUBLE PRECISION  A(N,N),B(N,N),TOL
      LOGICAL DONE
C
C   Local Declarations
       INTEGER*4 I,J
       DOUBLE PRECISION MNORM,RTOL
       LOGICAL MCRIT2
C
C  Preform A = A + b

       DO 20 J=1,N
         DO 10 I=1,J
          A(I,J) = A(I,J) + B(I,J)
          A(J,I) = A(I,J)
   10   CONTINUE
   20 CONTINUE


c  Check convergence 

      RTOL = TOL* MNORM(A,N,N)
      DONE = MCRIT2 (B,RTOL,N,N)



      RETURN
      END













