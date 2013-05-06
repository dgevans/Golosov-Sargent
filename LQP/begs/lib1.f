C_________________________________________________________________
C
C        LIBRARY:  FORTRAN LIBRARY I
C
C        SYNOPSIS:  Contains basic functions for copying matrices
C           and for checking the inputs to MEX functions 
C                
C        CONTENTS:
C            CHECKIN -   Checks the inputs to a MEX function for
C                         common errors
C            SCOPYI2R     Copies a scalar ineteger*4 to a scalar
C                          double precision variable
C            MCONSR  -   Brings a matrix stored in symmetric fromat to
C                         full storage format
C            MCOPYM  -   Copies a matrix
C            MCOPYM2 -   Copies a matrix when the row dimension
C                           of the resulting matrix is different
C                           from the matrix being copied
C            MCOPYS  -   Copies a symmteric matrix
C
C        WRITTEN BY:  Evan Anderson     January 1, 1995
C_________________________________________________________________



C*********************************************************************
C
C     SUBROUTINE:  CHECKIN
C
C     SYNOPSIS:  Checks the inputs to a MEX file to make sure none
C                of the inputs are matrices that are not
C                           (1) Of type double
C                           (2) Full
C                           (3) Numeric	
C
C*******************************************************************

      SUBROUTINE CHECKIN (NUMIN,PRHS)
      INTEGER*4  NUMIN,PRHS(*)
     
      INTEGER*4 NUMDOUB,NUMFULL,NUMNUM,J

      NUMDOUB =0
      NUMFULL = 0
      NUMNUM=0
      DO 10 J = 1,NUMIN
        NUMDOUB = NUMDOUB + MXISDOUBLE(PRHS(J))
        NUMFULL = NUMFULL + MXISFULL(PRHS(J))
        NUMNUM = NUMNUM + MXISNUMERIC (PRHS(J))
   10 CONTINUE
      IF (NUMDOUB .NE.NUMIN ) THEN
        CALL MEXERRMSGTXT('All inputs must be of type double')
      ELSEIF (NUMFULL .NE. NUMIN) THEN
        CALL MEXERRMSGTXT ('All inputs must be full matrices')
      ELSEIF (NUMNUM .NE. NUMIN) THEN
        CALL MEXERRMSGTXT ('All inputs must be numeric matrices')
      ENDIF

      RETURN
      END

C*********************************************************************
C
C     SUBROUTINE:  SCOPYI2R
C
C     SYNOPSIS:  Copies a scalar integer*4 variable to a double
C                 precision variable.  This is a very useful
C                 subroutine for passing integer variables out
C                 of MEX files
C
C*******************************************************************

      SUBROUTINE  SCOPYI2R(XINT,XREAL)
      DOUBLE PRECISION XREAL
      INTEGER*4 XINT
 
      XREAL = REAL(XINT)
 
      RETURN
      END


          


C**************************************************************
C
C     SUBROUTINE:  MCONSR
C
C     SYNOPSIS:  Brings a matrix stored in symmetric format to full
C                 storage format
C
C**************************************************************

      SUBROUTINE MCONSR(A,N)
      DOUBLE PRECISION A(N,N)
      INTEGER*4 N


      INTEGER I,J

      DO 20  I =1,N
       DO 10 J=1,I-1
          A(I,J) = A(J,I)
   10  CONTINUE
   20 CONTINUE
      RETURN
      END

C******************************************************************
C
C    SUBROUTINE:  MCOPYM
C
C    SYNOPSIS:  Copies a mxn matrix A into a mxn matrix B
C
C    PARAMETERS:
C      A       MxN double precision matrix
C      B       MxN double precision matrix
C      M       row dimension of A and B
C      N       column dimension of A and B 
C      
C
C******************************************************************

      SUBROUTINE MCOPYM (A,B,M,N)
      INTEGER*4 M,N
      DOUBLE PRECISION    A(M,N), B(M,N)
C
C Local Declarations      
      INTEGER*4  I,J
C
C Copy A into B using nested do loops
      DO 20 J=1,N
        DO 10 I=1,M
           B(I,J) = A(I,J)
   10   CONTINUE
   20 CONTINUE
C
C Finish
      RETURN
      END

C******************************************************************
C
C     SUBROUTINE:  MCOPYM2
C
C     SYNOPSIS:  Copies a mxn matrix A into B when the B is
C                  declared with a geater than m row
C                  dimension. elements outside out the mxn
C                  range in B are not altered
C
C*******************************************************************

      SUBROUTINE MCOPYM2(A,B,M,N,M2,N2)
      INTEGER*4 M,N,M2,N2
      REAL*8 A(M,N), B(M2,N2)
C
C This function copies the matrix A into the matrix B.
C It is assumed that M2.GE.M  and N2.GE.N.  This assumption
C    is not checked
C
C
C  Local variables
      INTEGER*4 I,J
C
C  The main loop
      DO 20 J=1,N
        DO 10 I=1,M
          B(I,J) = A(I,J)
   10  CONTINUE
   20 CONTINUE
C
C finish
      RETURN
      END



C******************************************************************
C
C    SUBROUTINE:  MCOPYS
C
C    SYNOPSIS:  Copies a nxn symmetric matrix A into a nxn symmetric matrix B
C       The fact that A is symmetric is not checked.  The strict lower
C       triangular portions of A and B are not refernced. Upon exit, the
C       upper triangular part of B will equal the upper triangular part of A.
C       The lower triangular portion of B will remain unchanged.
C
C
C    PARAMETERS:
C      A       NxN double precision matrix
C      B       NxN double precision matrix
C      N       order of A and B
C
C
C******************************************************************


      SUBROUTINE MCOPYS(A,B,N)
      INTEGER*4  N
      DOUBLE PRECISION  A(N,N), B(N,N)
C
C Local Parameters
      INTEGER*4   I,J
C
C Copy the upper triangular portion of A into B
      DO 20 J=1,N
        DO 10 I=1,J    
          B(I,J) = A(I,J) 
   10   CONTINUE
   20 CONTINUE
C
C Finish
      RETURN
      END

























    
