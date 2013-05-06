C_________________________________________________________________
C
C        LIBRARY:  FORTRAN LIBRARY 4
C
C        SYNOPSIS:  Contains functions for Hessenberg matrices
C
C        CONTENTS:
C          DRIVER ROUTINES
C            FULL2HESS-  Copies a full matrix to a matrix in Hessenberg
C                          storage
C            HESS2FULL - Inverse of FULL2HESS
C            VLINEARH-   Solves the linear system Ax=b when A is hessenberg 
C          FROM ALGORITHM 705 of the ACM
C            HSFA
C            HSSL
C
C_________________________________________________________________



C******************************************************************************
C
C     SUBROUTINE:  FULL2HESS
C
C     SYNOPSIS:   Copies the NxN full matrix F   into a vector
C                 Hessenberg representation.
C
C                      
C
C*******************************************************************************



      SUBROUTINE FULL2HESS(F,H,N)
      DOUBLE PRECISION F(N,N), H(*)
      INTEGER*4 N

      INTEGER*4 I,J,K

c  handle all of F execpt for the last column
      K=0
      DO 20 J=1,N-1
        DO 10 I = 1,J+1
          K=K+1
          H(K) = F(I,J)
   10   CONTINUE
   20 CONTINUE    
c
c handle the last column
      J = N
      DO 30 I=1,N
         K = K +1
         H(K) = F(I,J)
   30 CONTINUE
c
c  finish
      RETURN
      END



C******************************************************************************
C
C     SUBROUTINE:  HESS2FULL 
C
C     SYNOPSIS:    Copies a matrix H in vectot Hessenberg storage into
C                    a full matrix F
C
C*******************************************************************************




      SUBROUTINE HESS2FULL(H,F,N)
      DOUBLE PRECISION F(N,N), H(*)
      INTEGER*4 N

      INTEGER*4 I,J,K

c
c  handle all of F execpt for the last column
      K=0
      DO 20 J=1,N-1
        DO 10 I = 1,J+1
          K=K+1
          F(I,J) = H(K)
   10   CONTINUE
   20 CONTINUE
c
c handle the last column
      J = N
      DO 30 I=1,N
         K = K +1
         F(I,J) = H(K)
   30 CONTINUE
c
c  finish
      RETURN
      END



C***************************************************************
C
C     SUBROUTINE : VLINEARH
C
C     SYNOPSIS :  Solves two linear systems, Ax=b and
C                   when b is a vector.  A is assumed 
C                   to be a Hessenberg matrix
C
C     PARAMETERS
C          A     ON INPUT:  MxM double precision matrix
C                ON OUTPUT: lu decompsoition of A
C          B     ON INPUT: Mx1
C                ON OUTPUT  the solution x
C
C***************************************************************

      SUBROUTINE VLINEARH(A,B,M,WORKI,ERR)
      INTEGER*4 M, WORKI(M)
      DOUBLE PRECISION  A(M,M), B(M)
      LOGICAL ERR

c
c  local declarations
      INTEGER*4  INFO
c
c  get the lu decompsoition
      CALL HSFA (A,M,WORKI,INFO,1)
      ERR = (INFO.NE.0)
      IF (ERR) RETURN
c
c  solve Ax=b
      CALL HSSL(A,M,WORKI,B,1)
c
c  finish
   10 CONTINUE
      RETURN
      END       

c THIS IS A MODIFIED ROUTINE DESIGNED FOR 
c USE WITH MATLAB MEX FILES - THE ONLY CHANGE IS THAT INTEGER
c VARIABLES WERE CHANGED TO INTEGER*4 VARIABLES
c PROGRAM MODIFIED Nov 1994 by Evan Anderson
      SUBROUTINE HSFA(AV,N,IPVT,INFO,SD)
C
      INTEGER*4 N,IPVT(1),INFO,SD
      DOUBLE PRECISION AV(1)
C
C     THIS ROUTINE FACTORS A DOUBLE PRECISION MATRIX BY GAUSSIAN
C     ELIMINATION.  IT IS A MODIFIED VERSION OF DGEFA IN WHICH THE
C     MATRIX A HAS THE STRUCTURE OF AN UPPER TRIANGULAR MATRIX PLUS
C     SD NON-ZERO SUBDIAGONALS.  TO SAVE SPACE THE TWO-DIMENSIONAL
C     MATRIX A IS REPRESENTED AS A ONE-DIMENSIONAL ARRAY AV.  EACH
C     SUCCESSIVE COLUMN OF A IS STORED IN AV, IN SUCH A WAY THAT THE
C     J'TH COLUMN OF A IS ALLOTTED J+SD POSITIONS IN AV .
C     NOTE: IN THE CALLING PROGRAM THE DIMENSION OF IPVT AND Z SHOULD BE
C           AT LEAST N .  SET THE DIMENSION OF AV TO AT LEAST
C           (N*(N+1))/2 + N*SD .
C
C     REFS: J.J. DONGARRA, J.R. BUNCH, C.B. MOLER, AND G.W. STEWART,
C           LINPACK USERS' GUIDE, SIAM, 1979.
C
C     ON ENTRY -
C        AV      DOUBLE PRECISION (N*(N+1)/2 + N*SD)
C                ONE-DIMENSIONAL ARRAY REPRESENTING THE MATRIX A TO BE
C                FACTORED.
C
C        N       INTEGER*4
C                THE ORDER OF THE MATRIX  A .
C
C        SD      INTEGER*4
C                THE NUMBER OF NON-ZERO SUBDIAGONALS OF A .
C
C     ON RETURN -
C        AV      DOUBLE PRECISION (N*(N+1)/2)
C                A VECTOR REPRESENTING AN UPPER TRIANGULAR MATRIX AND
C                THE MULTIPLIERS USED TO OBTAIN IT.  THE FACTORIZATION
C                CAN BE WRITTEN  A = L*U  WHERE  L  IS A PRODUCT OF
C                PERMUTATION AND UNIT LOWER TRIANGULAR MATRICES AND
C                U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER*4(N)
C                AN INTEGER*4 VECTOR OF PIVOT INDICES.
C
C        INFO    INTEGER*4
C                = 0  NORMAL VALUE.
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                     INDICATE THAT HSSL WILL DIVIDE BY ZERO IF
C                     CALLED.  USE RCOND IN HSCO FOR A RELIABLE
C                     INDICATION OF SINGULARITY.
C
C     MODIFIED VERSION OF LINPACK ROUTINE DGEFA,  J. AMATO, APRIL 1984.
C     REVISED -
C       28NOV88 J. GARDINER, OSU CIS, COLUMBUS, OH  43210  (614)292-8658
c THIS IS A MODIFIED ROUTINE DESIGNED FOR 
c USE WITH MATLAB MEX FILES - THE ONLY CHANGE IS THAT INTEGER
c VARIABLES WERE CHANGED TO INTEGER*4 VARIABLES
c PROGRAM MODIFIED Nov 1994 by Evan Anderson
C
C     SUBROUTINES AND FUNCTIONS CALLED -
C       (LINPACK) DAXPY DSCAL IDAMAX
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION T
      INTEGER*4 IDAMAX,J,K,KP1,L,NM1,I1,I2,I3
C
      IF (N-SD .LT. 1) RETURN
C
C     INSERT "DUMMY" ELEMENTS FOR EASE OF SUBSEQUENT CODING
C
      DO 8 J = N-SD+1,N
         I1 = SD*(J-1) + (J*(J-1))/2
         DO 5 I2 = N+1,J+SD
            AV(I1+I2) = 0.0D0
    5    CONTINUE
    8 CONTINUE
C
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
*
C
C        FIND L = PIVOT INDEX
C
         I1 = K + SD*(K-1) + (K*(K-1))/2
         L = IDAMAX(SD+1,AV(I1),1) + K - 1
         IPVT(K) = L
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         I2 = I1 + L - K
         IF (AV(I2) .EQ. 0.0D0) GO TO 40
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .EQ. K) GO TO 10
               T = AV(I2)
               AV(I2) = AV(I1)
               AV(I1) = T
   10       CONTINUE
C
C           COMPUTE MULTIPLIERS
C
            T = -1.0D0/AV(I1)
            CALL DSCAL(SD,T,AV(I1+1),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            DO 30 J = KP1, N
               I2 = L + SD*(J-1) + (J*(J-1))/2
               I3 = I2 + K - L
               T = AV(I2)
               IF (L .EQ. K) GO TO 20
                  AV(I2) = AV(I3)
                  AV(I3) = T
   20          CONTINUE
               CALL DAXPY(SD,T,AV(I1+1),1,AV(I3+1),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      I1 = N + SD*(N-1) + (N*(N-1))/2
      IF (AV(I1) .EQ. 0.0D0) INFO = N
      RETURN
C --- LAST LINE OF HSFA ---
      END


c THIS IS A MODIFIED ROUTINE DESIGNED FOR 
c USE WITH MATLAB MEX FILES - THE ONLY CHANGE IS THAT INTEGER
c VARIABLES WERE CHANGED TO INTEGER*4 VARIABLES
c PROGRAM MODIFIED Nov 1994 by Evan Anderson
      SUBROUTINE HSSL(AV,N,IPVT,B,SD)
C
      INTEGER*4 N,IPVT(1),SD
      DOUBLE PRECISION AV(1),B(1)
C
C     THIS ROUTINE SOLVES THE DOUBLE PRECISION SYSTEM
C     A * X = B  OR  TRANS(A) * X = B
C     USING THE FACTORS COMPUTED BY HSCO OR HSFA.
C     THE INPUT MATRIX IS IN THE FORM OF AN UPPER TRIANGULAR MATRIX
C     PLUS SD NON-ZERO SUBDIAGONALS.  TO SAVE SPACE THE TWO-DIMENSIONAL
C     MATRIX A IS REPRESENTED AS A ONE-DIMENSIONAL ARRAY AV.  EACH
C     SUCCESSIVE COLUMN OF A IS STORED IN AV, SUCH THAT THE J'TH COLUMN
C     OF A IS ALLOTTED J+SD POSITIONS IN AV .
C     NOTE: IN THE CALLING PROGRAM THE DIMENSION OF IPVT AND B SHOULD BE
C           AT LEAST N .  SET THE DIMENSION OF AV TO AT LEAST
C           (N*(N+1))/2 + N*SD .
C
C     REFS: J.J. DONGARRA, J.R. BUNCH, C.B. MOLER, AND G.W. STEWART,
C           LINPACK USERS' GUIDE, SIAM, 1979.
C
C     ON ENTRY -
C
C        AV      DOUBLE PRECISION (N*(N+1)/2 + SD)
C                ONE-DIMENSIONAL ARRAY REPRESENTING THE MATRIX A (OUTPUT
C                FROM HSCO OR HSFA).
C
C        N       INTEGER*4
C                THE ORDER OF THE MATRIX  A .
C
C        B       DOUBLE PRECISION VECTOR
C                DATA VECTOR (RIGHT HAND SIDE).
C
C        SD      INTEGER*4
C                THE NUMBER OF NON-ZERO SUBDIAGONALS OF A .
C
C        IPVT    INTEGER*4(N)
C                THE PIVOT VECTOR FROM HSCO OR HSFA.
C
C     ON RETURN -
C
C        B       CONTAINS THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION -
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
C        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS.  IT WILL NOT
C        OCCUR IF THE SUBROUTINES ARE CALLED CORRECTLY AND IF HSCO HAS
C        SET RCOND .GT. 0.0 OR HSFA HAS SET INFO .EQ. 0 .
C
C     MODIFIED VERSION OF LINPACK ROUTINE DGESL, J. AMATO, APRIL 1984.
C     REVISED -
C       28NOV88 J. GARDINER, OSU CIS, COLUMBUS, OH  43210  (614)292-8658
c THIS IS A MODIFIED ROUTINE DESIGNED FOR 
c USE WITH MATLAB MEX FILES - THE ONLY CHANGE IS THAT INTEGER
c VARIABLES WERE CHANGED TO INTEGER*4 VARIABLES
c PROGRAM MODIFIED Nov 1994 by Evan Anderson
C
C     SUBROUTINES AND FUNCTIONS CALLED -
C       (LINPACK) DAXPY DDOT
C
C     INTERNAL VARIABLES -
C
      DOUBLE PRECISION T
      INTEGER*4 J,K,KB,L,NM1,I1
C
      IF (N-SD .LT. 1) RETURN
C
C     INSERT "DUMMY" ELEMENTS FOR EASE OF SUBSEQUENT CODING
C
      DO 8 J = N-SD+1,N
         I1 = SD*(J-1) + (J*(J-1))/2
         DO 5 L = N+1,J+SD
            AV(I1+L) = 0.0D0
    5    CONTINUE
    8 CONTINUE
C
      NM1 = N - 1
C
C        FIRST SOLVE  L*Y = B
C
         IF (NM1 .LT. 1) GO TO 30
         DO 20 K = 1, NM1
            L = IPVT(K)
            T = B(L)
            IF (L .EQ. K) GO TO 10
               B(L) = B(K)
               B(K) = T
   10       CONTINUE
            I1 = K + SD*(K-1) + (K*(K-1))/2
            CALL DAXPY(SD,T,AV(I1+1),1,B(K+1),1)
   20    CONTINUE
   30    CONTINUE
C
C        NOW SOLVE  U*X = Y
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            I1 = K + SD*(K-1) + (K*(K-1))/2
            B(K) = B(K)/AV(I1)
            T = -B(K)
            CALL DAXPY(K-1,T,AV(I1+1-K),1,B(1),1)
   40    CONTINUE
      RETURN
C --- LAST LINE OF HSSL ---
      END
