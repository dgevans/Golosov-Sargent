C-------------------------------------------------------------------------------
C        LIBRARY:  FORTRAN LIBRARY 2: Matrix Arithmetic
C
C        SYNOPSIS:  Contains functions for arithmetic operations on
C           symmetric and non-symmetric matrices
C
C        CONTENTS:
C            MADD      -  A = A + B
C            MADD3     -  A = A + B + C
C            MCRIT     -  returns true if every column sum of a 
C                           of a matrix  is less than some tolerance
C            MCRIT2    -  returns true if every column sum of a matrix 
C                         in symmetric storage is less than some tolerance
C            MCRITR    -  a version of MCRIT for relative tolerance
C            MINCI     -  A = A + I   where I is the indentity matrix
C            MMUL      -  C = A*B  
C            MMUL2     -  C = A'*B
C            MMUL3     -  C = A*B'
C            MMULADDW  -  C = C + A*B  
C            MMULADDW2 -  D = C + A*B       
C            MMULMMS   -  C = A*B   where C is symmetric.  Only upper tringular part
C                           of C is calculated 
C            MMULMMSF  -   C = A*B   where C is  symmetric.  All of C is calculated
C            MMULQUAD  -  C = A'*B*A where C is symmetric 
C            MMULSM    -  C = A*B    where A is symmetric
C            MMULSUBW2 -  D = C - A*B
C            MNORM     -  Returns the maximum column  sum of a matrix
C            MNORM2    -   Returns the meximum coulmn sum for a matrix
C                           in symmetric storage
C
C        WRITTEN BY:   Evan Anderson     November 15, 1994
C        MODIFIED:    12/27/98 
C_________________________________________________________________


C*********************************************************************
C
C    SUBROUTINE: MADD
C
C    SYNOPSIS:  Performs the matrix addition A = A + B
C
C    PARAMETERS:
C       A   MxN double precision matrix
C       B   MxN double precision matrix
C       M   row dimension of A and B
C       N   column dimension of A and B
C
C*********************************************************************

      SUBROUTINE MADD(A,B,M,N)
      INTEGER*4  M,N
      DOUBLE PRECISION A(M,N),B(M,N)
C
C  Local Declarations
      INTEGER*4 I,J
C
C  The main loop
      DO 20 J=1,N
        DO 10 I=1,M
          A(I,J) = A(I,J) + B(I,J)
   10   CONTINUE
   20 CONTINUE
C
C finish
      RETURN
      END


C*********************************************************************
C
C    SUBROUTINE: MADD3
C
C    SYNOPSIS:  Performs the matrix addition  C =  A + B + C
C       where {A,B,C}  are M x N matrices
C
C*********************************************************************

      SUBROUTINE MADD3(A,B,C,M,N)
      INTEGER*4  M,N
      DOUBLE PRECISION A(M,N),B(M,N),C(M,N)
C
C  Local Declarations
      INTEGER*4 I,J
C
C  The main loop
      DO 20 J=1,N
        DO 10 I=1,M
          C(I,J) = A(I,J) + B(I,J) + C (I,J)
   10   CONTINUE
   20 CONTINUE
C
C finish
      RETURN
      END 




C*******************************************************************
C
C   FUNCTION:  MCRIT   
C
C   SYNOPSIS:  Returns true if every column sum of the MxN matrix
C       A is less than  tol
C
C********************************************************************

      LOGICAL FUNCTION  MCRIT (A,TOL,M,N)
      DOUBLE PRECISION A(M,N), TOL
      INTEGER*4  M,N


      DOUBLE PRECISION SUM,DABS
      INTEGER*4  I,J

      DO 20 J=1,N
        SUM = 0
        DO 10 I=1,M 
          SUM = SUM + DABS(A(I,J)) 
   10   CONTINUE
        IF (SUM .GT. TOL) THEN
           MCRIT = .FALSE.
           RETURN
        ENDIF
   20 CONTINUE
      MCRIT = .TRUE.
      RETURN
      END


C*******************************************************************
C
C   FUNCTION:  MCRIT2 
C
C   SYNOPSIS:  A version of MCRIT when A is symmetric 
C
C********************************************************************
 
      LOGICAL FUNCTION  MCRIT2 (A,TOL,M,N)
      DOUBLE PRECISION A(M,N), TOL
      INTEGER*4  M,N
 
 
      DOUBLE PRECISION SUM,DABS
      INTEGER*4  I,J
 
      DO 20 J=1,N
        SUM = 0
        DO 10 I=1,J
          SUM = SUM + DABS(A(I,J))
   10   CONTINUE
        DO 15 I=J+1,M
          SUM = SUM + DABS(A(J,I))
   15   CONTINUE
        IF (SUM .GT. TOL) THEN
           MCRIT2 = .FALSE.
           RETURN
        ENDIF
   20 CONTINUE
      MCRIT2 = .TRUE.
      RETURN
      END



C*******************************************************************
C
C   FUNCTION:  MCRITR
C
C   SYNOPSIS:  Returns true if every column sum of A is < norm(B)*tol
C
C********************************************************************

      LOGICAL FUNCTION  MCRITR (A,B,TOL,M,N)
      DOUBLE PRECISION A(M,N), B(M,N),TOL
      INTEGER*4  M,N

      LOGICAL MCRIT
      DOUBLE PRECISION MNORM

      
      MCRITR = MCRIT (A,TOL* MNORM(B,M,N),M,N)
      RETURN
      END


C*******************************************************************
C
C    SUBROUTINE: MINCI
C
C    SYNOPSIS:  Perfoms the matrix addition A = A + Identity_matrix
C               A is overwritten
C
C    PARAMATERS
C        A   NxN double precision matrix
C        N   column and row dimension of A
C
C*********************************************************************      

      SUBROUTINE MINCI(A,N)
      INTEGER*4  N
      DOUBLE PRECISION  A(N,N)

C   Local Declarations
      INTEGER*4 I
C
C  The main loop
      DO 10 I=1,N
        A(I,I) = A(I,I) + 1
   10 CONTINUE
C
C finish
      RETURN
      END





C******************************************************************
C
C   SUBROUTINE: MMUL
C
C   SYNOPSIS: Computes the matrix product C = A*B
C
C   PARAMETERS:
C     A   MxN double precision matrix
C     B   N*P double precision matrix
C     C   MxP double precision matrix
C     M   row dimension of A and C
C     N   column dimension of A and row dimension of B
C     P   column dimension of B and column dimension of P
C
C*******************************************************************


      SUBROUTINE MMUL(A,B,C,M,N,P)
      INTEGER*4  M,N,P
      DOUBLE PRECISION  A(M,N),B(N,P),C(M,P)
C
C   Local Declarations
      INTEGER*4 I,J,K
      DOUBLE PRECISION  TEMP
C
C  The main loop
      DO 30 I=1,M
        DO 20 J=1,P
          TEMP =0
          DO 10 K=1,N
            TEMP = TEMP+A(I,K)*B(K,J)
   10     CONTINUE
          C(I,J) = TEMP
   20   CONTINUE
   30 CONTINUE
C
C finish
      RETURN
      END


C******************************************************************
C
C   SUBROUTINE: MMUL2
C
C   SYNOPSIS: Computes the matrix product C = A'*B
C
C   PARAMETERS:
C     A   N*M double precision matrix
C     B   N*P double precision matrix
C     C   MxP double precision matrix
C     M   column dimension of A and row dimension of C
C     N   row dimension of A and row dimension of B
C     P   column dimension of B and coulmn dimesnion of C
C
C*******************************************************************
 
 
      SUBROUTINE MMUL2(A,B,C,M,N,P)
      INTEGER*4  M,N,P
      DOUBLE PRECISION  A(N,M),B(N,P),C(M,P)
C
C   Local Declarations
      INTEGER*4 I,J,K
      DOUBLE PRECISION  TEMP
C
C  The main loop
      DO 30 I=1,M
        DO 20 J=1,P
          TEMP =0
          DO 10 K=1,N
            TEMP = TEMP+A(K,I)*B(K,J)
   10     CONTINUE
          C(I,J) = TEMP
   20   CONTINUE
   30 CONTINUE
C
C finish
      RETURN
      END

C******************************************************************
C
C   SUBROUTINE: MMUL3
C
C   SYNOPSIS: Computes the matrix product C = A*B'
C
C   PARAMETERS:
C     A   MxN double precision matrix
C     B   PxN double precision matrix
C     C   MxP double precision matrix
C     M   row dimension of C and A
C     N   column dimension of C and A
C     P   row dimension of B and column dimension of C
C
C*******************************************************************
 
 
      SUBROUTINE MMUL3(A,B,C,M,N,P)
      INTEGER*4  M,N,P
      DOUBLE PRECISION  A(M,N),B(P,N),C(M,P)
C
C   Local Declarations
      INTEGER*4 I,J,K
      DOUBLE PRECISION  TEMP
C
C  The main loop
      DO 30 I=1,M
        DO 20 J=1,P
          TEMP =0
          DO 10 K=1,N
            TEMP = TEMP+A(I,K)*B(J,K)
   10     CONTINUE
          C(I,J) = TEMP
   20   CONTINUE
   30 CONTINUE
C
C finish
      RETURN
      END

C*****************************************************************************
C
C     SUBROUTINE:  MMULADDW
C
C     SYNOPSIS:  Performs the operation C = C + A*B where A is MxN and B is NxP
C
C*******************************************************************************



      SUBROUTINE MMULADDW(A,B,C,M,N,P,LA,LB,LC)
      INTEGER*4  M,N,P,LA,LB,LC
      DOUBLE PRECISION  A(LA,N),B(LB,P),C(LC,P)
C
C   Local Declarations
      INTEGER*4 I,J,K
      DOUBLE PRECISION  TEMP
C
C  The main loop
      DO 30 I=1,M
        DO 20 J=1,P
          TEMP =0
          DO 10 K=1,N
            TEMP = TEMP+A(I,K)*B(K,J)
   10     CONTINUE
          C(I,J) = C(I,J)+TEMP
   20   CONTINUE
   30 CONTINUE
C
C finish
      RETURN
      END


C*****************************************************************************
C
C     SUBROUTINE:  MMULADDW2
C
C     SYNOPSIS:  Performs the operation D = C + A*B where A is MxN and B is NxP
C
C*******************************************************************************

      SUBROUTINE MMULADDW2(A,B,C,D,M,N,P,LA,LB,LC,LD)
      INTEGER*4  M,N,P,LA,LB,LC,LD
      DOUBLE PRECISION  A(LA,N),B(LB,P),C(LC,P),D(LD,P)
C
C   Local Declarations
      INTEGER*4 I,J,K
      DOUBLE PRECISION  TEMP
C
C  The main loop
      DO 30 I=1,M
        DO 20 J=1,P
          TEMP =0
          DO 10 K=1,N
            TEMP = TEMP+A(I,K)*B(K,J)
   10     CONTINUE
          D(I,J) = C(I,J)+TEMP
   20   CONTINUE
   30 CONTINUE
C
C finish
      RETURN
      END


C*****************************************************************************
C
C     SUBROUTINE:  MMULSUBW2
C
C     SYNOPSIS:  Performs the operation D = C - A*B where A is MxN and B is NxP
C
C******************************************************************************* 
      SUBROUTINE MMULSUBW2(A,B,C,D,M,N,P,LA,LB,LC,LD)
      INTEGER*4  M,N,P,LA,LB,LC,LD
      DOUBLE PRECISION  A(LA,N),B(LB,P),C(LC,P),D(LD,P)
C
C   Local Declarations
      INTEGER*4 I,J,K
      DOUBLE PRECISION  TEMP
C
C  The main loop
      DO 30 I=1,M
        DO 20 J=1,P
          TEMP =0
          DO 10 K=1,N
            TEMP = TEMP+A(I,K)*B(K,J)
   10     CONTINUE
          D(I,J) = C(I,J)-TEMP
   20   CONTINUE
   30 CONTINUE
C
C finish
      RETURN
      END
 



C******************************************************************
C
C    SUBROUTINE:  MMULMMS
C
C    SYNOPSIS:  Computes the matrix product C = AxB , when C is 
C        assumed to be symmetric. Only the upper triangular portion
C        of C = AxB is calculated. C is the only parameter that is
C        changed upon exit 
C
C    PARAMETERS:
C      A    MxN double precision matrix
C      B    NxM double precision matrix
C      C    MxM double precision matrix
C      M    Row dimension of A and C.  Column dimension of B and C.
C      N    Column dimension of A.  Row dimension of B. 
C
C******************************************************************


      SUBROUTINE MMULMMS(A,B,C,M,N)
      INTEGER*4  M,N
      DOUBLE PRECISION  A(M,N),B(N,M),C(M,M)
C
C   Local Declarations
      INTEGER*4 I,J,K
      DOUBLE PRECISION  TEMP
C
C  The main loop
      DO 30 I=1,M
        DO 20 J=I,M
          TEMP =0
          DO 10 K=1,N
            TEMP = TEMP+A(I,K)*B(K,J)
   10     CONTINUE
          C(I,J) = TEMP
   20   CONTINUE
   30 CONTINUE
C
C finish
      RETURN
      END



C******************************************************************
C
C    SUBROUTINE:  MMULQUAD
C
C    SYNOPSIS:  Computes the matrix product C = A'*B*A where
C        symmetric storage is used for the symmetric matrix C.   
C
C    PARAMETERS:
C      A    MxN double precision matrix
C      B    MxM double precision matrix
C      C    NxN double precision matrix
C      WORK  NxM double precision matrix
C
C******************************************************************


      SUBROUTINE MMULQUAD(A,B,C,WORK,M,N)
      INTEGER*4  M,N
      DOUBLE PRECISION  A(M,N),B(M,M),C(N,N),WORK(N,M)
C
c
c  do the multiplications
      CALL MMUL2 (A,B,WORK,N,M,M)
      CALL MMULMMS (WORK,A, C, N,M)
c  
c  finish
      RETURN
      END






C******************************************************************
C
C    SUBROUTINE:  MMULMMSF
C
C    SYNOPSIS: Same as MMULMMS.  C is symmetric.  But here we
C        store C in full format--- both the upper and lower
C        triangular portions of C are changed on exit  
C
C    PARAMETERS:
C      A    MxN double precision matrix
C      B    NxM double precision matrix
C      C    MxM double precision matrix
C      M    Row dimension of A and C.  Column dimension of B and C.
C      N    Column dimension of A.  Row dimension of B.
C
C******************************************************************


      SUBROUTINE MMULMMSF(A,B,C,M,N)
      INTEGER*4  M,N
      DOUBLE PRECISION  A(M,N),B(N,M),C(M,M)
C
C   Local Declarations
      INTEGER*4 I,J,K
      DOUBLE PRECISION  TEMP
C
C  The main loop
      DO 30 I=1,M
        DO 20 J=I,M
          TEMP =0
          DO 10 K=1,N
            TEMP = TEMP+A(I,K)*B(K,J)
   10     CONTINUE
          C(I,J) = TEMP
          C(J,I) = TEMP
   20   CONTINUE
   30 CONTINUE
C
C finish
      RETURN
      END


C******************************************************************
C
C    SUBROUTINE:  MMULSM
C
C    SYNOPSIS:  Computes the matrix product C = AxB , when A is
C        assumed to be symmetric. Only the upper triangular portion
C        of A is referenced. C is the only parameter that is
C        changed upon exit
C
C    PARAMETERS:
C      A    NxN double precision matrix
C      B    NxP double precision matrix
C      C    NxP double precision matrix
C      N    Order of A.  Row dimension of B and C
C      P    Column dimension of B and C
C
C********************************************************************

      SUBROUTINE MMULSM(A,B,C,N,P)
      INTEGER*4  N,P
      DOUBLE PRECISION  A(N,N),B(N,P),C(N,P)
C      C    NxP double precision matrix
C      N    order of A;  Number of rows in B and C
C      P    number of columns in B and C 
C
C   Local Declarations
       INTEGER*4 I,J,K
       DOUBLE PRECISION  TEMP
C
C  The main loop
      DO 30 I=1,N
        DO 20 J=1,P
          TEMP =0
          DO 5 K=1,I
            TEMP = TEMP + A(K,I)*B(K,J)
    5     CONTINUE      
          DO 10 K=I+1,N
            TEMP = TEMP+A(I,K)*B(K,J)
   10     CONTINUE
          C(I,J) = TEMP
   20   CONTINUE
   30 CONTINUE
C
C finish
      RETURN
      END



C*******************************************************************
C
C   FUNCTION:  MNORM    
C
C   SYNOPSIS:  Returns the maximum  (in absolute value) column
C            sum of the matrix A 
C
C********************************************************************
 
      DOUBLE PRECISION  FUNCTION MNORM (A,M,N)
      DOUBLE PRECISION A(M,N)
      INTEGER*4  M,N
 
 
      DOUBLE PRECISION SUM,DABS
      DOUBLE PRECISION BEST
      INTEGER*4  I,J

      BEST = 0
      DO 20 J=1,N
        SUM = 0
        DO 10 I=1,M
          SUM = SUM + DABS(A(I,J))
   10   CONTINUE
        IF (SUM .GT. BEST) BEST = SUM
   20 CONTINUE
      MNORM = BEST
      RETURN
      END





C*******************************************************************
C
C   FUNCTION:  MNORM2    
C
C   SYNOPSIS:   Returns the maximum  (in absolute value) column
C               sum of the matrix A when A is symmetric
C
C
C********************************************************************

      DOUBLE PRECISION FUNCTION MNORM2 (A,M,N)
      DOUBLE PRECISION A(M,N)
      INTEGER*4  M,N

      DOUBLE PRECISION  SUM,DABS
      DOUBLE PRECISION BEST
      INTEGER*4  I,J

      BEST = 0
      DO 20 J=1,N
        SUM = 0
        DO 10 I = 1,J
           SUM= SUM + DABS(A(I,J))
   10   CONTINUE
        DO 15  I = J+1,M
           SUM = SUM + DABS(A(J,I))
   15   CONTINUE
      IF (SUM .GT. BEST) BEST = SUM
   20 CONTINUE
      MNORM2 = BEST
      RETURN
      END 







