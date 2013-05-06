C_________________________________________________________________
C
C        LIBRARY:  FORTRAN LIBRARY 3:  LINPACK ROUTINES
C
C        SYNOPSIS: driver programs for subroutines from
C           LINPACK 
C
C        CONTENTS:
C           (1) Driver routines
C                  MDET
C                  MDETINVS
C                  MLINEAR
C                  MLINEAR
C                  MLINEAR2
C                  VLINEAR
C           (2) LINPACK routines
C                  DGEDI
C                  DGEFA
C                  DGESL
C                  DSIDI
C                  DSIFA  
C
C__________________________________________________________________




C__________________________________________________________________
C
C      DRIVER ROUTINES
C
C__________________________________________________________________





C****************************************************************
C
C   SUBROUTINE :  MDET
C
C****************************************************************      


      SUBROUTINE MDET(A,N,DETERM,WORK,WORKR,WORKI,ERR)
      DOUBLE PRECISION A(N,N), DETERM
      DOUBLE PRECISION WORK(N,N),WORKR(N)   
      INTEGER*4  N      
      INTEGER*4 WORKI(N)
      LOGICAL ERR
C
C  Calculates the determinant of the NxN matrix A
C    without overwriting A
C
C  INPUTS 
C
C     A         NxN double precision matrix whose determinant
C                  is to be calculated
C     N         Order of A
C     DETERM    double precision scalar
C     WORK      NxN double precision work matrix
C     WORKR     N dimensional double precision work vector
C     WORKI     N dimensional double precision integer vector
C     ERR       Logical variable
C
C   OUTPUTS 
C
C     DETERM   The determinant of A
C     ERR      On normal return is false.  If A is singular, then
C                the determinant could not be calculated and ERR
C                is false
C
C
C
       DOUBLE PRECISION DET(2)
       INTEGER*4 INFO

C
C  Copy the input A and take the LU decomposition of A
       CALL MCOPYM (A,WORK,N,N)
       CALL DGEFA(WORK,N,N,WORKI,INFO)
       ERR = (INFO.NE.0)
       IF (ERR) GO TO 100
C
C  Calculate the determinant
       CALL DGEDI(WORK,N,N,WORKI,DET,WORKR,10)
       DETERM = DET(1) * 10.0**DET(2)        
C
C finish
  100 RETURN
      END
C
C
C*****************************************************************
C
C   SUBROUTINE : MDETINVS
C
C*****************************************************************

      SUBROUTINE MDETINVS(A,N,DETERM,WORKR,WORKI,ERR)
      DOUBLE PRECISION A(N,N), DETERM
      DOUBLE PRECISION WORKR(N)
      INTEGER*4  N
      INTEGER*4 WORKI(N)
      LOGICAL ERR
C
C  Calculates the inverse and determinant of the NxN symmetric
C     matrix A.  The fact that A is symmetric is not checked. The strict
C     lower trinagular portion of A is not referenced.  
C     A is overwriiten with its inverse. 
C
C  INPUTS
C
C     A         NxN double precision matrix whose determinant
C                  is to be calculated
C     N         Order of A
C     DETERM    double precision scalar
C     WORKR     N dimensional double precision work vector
C     WORKI     N dimensional double precision integer vector
C     ERR       Logical variable
C
C   Outputs
C
C     DETERM   The determinant of A
C     A         The inverse of A
C     ERR      On normal return is false.  If A is singular, then
C                the determinant could not be calculated and ERR
C                is false
C
       DOUBLE PRECISION DET(2)
       INTEGER*4 INFO, INERT(3)

C
C  Copy the input A and take the LU decomposition of A
       CALL DSIFA(A,N,N,WORKI,INFO)
       ERR = (INFO.NE.0)
       IF (ERR) GO TO 100
C
C  Calculate the determinant
       CALL DSIDI(A,N,N,WORKI,DET,INERT,WORKR,011)
       DETERM = DET(1) * 10**DET(2)
C
C finish
  100 RETURN
      END


C*********************************************************************
C
C    SUBROUTINE:  MINV
C
C    SYNOPSIS:  Calculates the inverse of a general matrix.  IF A
C        is singular then the inverse is not calculated and ERR=.TRUE..
C        A is overwritten with its inverse
C
C    PARAMETERS:
C        A       ON INPUT:   NxN double precision matrix
C                ON OUTPUT:  Inverse of the input
C        N       OREDER of A
C        WORKI   Integer working vector of length N
C        WORKR   Double precision working vector of length N
C        ERR     Boolean variable indicating failure
C
C***********************************************************************

      SUBROUTINE MINV(A,N,WORKI,WORKR,ERR)
      INTEGER*4  N,WORKI(N)
      DOUBLE PRECISION  A(N,N),WORKR(N)
      LOGICAL ERR


C   Local Declarations
       INTEGER*4 INFO
       DOUBLE PRECISION DET(2)
C
C  Triangularize
       CALL DGEFA(A,N,N,WORKI,INFO)
       ERR = (INFO.NE.0)
C
C  solve
       IF (.NOT. ERR) THEN
         CALL DGEDI(A,N,N,WORKI,DET,WORKR,01)
       ENDIF
C
C finish
      RETURN
      END

C***********************************************************
C
C    SUBROUTINE: MLINEAR  
C
C    SYNOPSIS:  Given A in its lu form, MLINAER solves the 
C       matrix linear equation AX=B.  Note that A must not
C       be singular.  A is MxM and B is MxN
C
C***********************************************************

      SUBROUTINE  MLINEAR(A,B,M,N,PIVOT)
      INTEGER*4   M,N, PIVOT(M)
      DOUBLE PRECISION A(M,M), B(M,N)

C
C Local variables
      INTEGER*4 J 

      DO 10 J = 1,N
         CALL DGESL(A,M,M,PIVOT,B(1,J),0)
   10 CONTINUE
      RETURN
      END



C***************************************************************
C
C     SUBROUTINE : MLINEAR1
C
C     SYNOPSIS :  Solves one linear system, AX=B 
C
C     PARAMETERS
C          A     ON INPUT:  MxM double precision matrix
C                ON OUTPUT: lu decompsoition of A
C          B     ON INPUT: MxN
C                ON OUTPUT  the solution X
C          C     ON INPUT:  MxN
C                ON OUTPUT  the solution Z
C
C***************************************************************

      SUBROUTINE MLINEAR1(A,B,M,N,WORKI,ERR)
      INTEGER*4 M,N, WORKI(M)
      DOUBLE PRECISION  A(M,M), B(M,N)
      LOGICAL ERR

c
c  local declarations
      INTEGER*4  INFO
c
c  get the lu decompsoition
      CALL DGEFA (A,M,M,WORKI,INFO)
      ERR = (INFO.NE.0)
      IF (ERR) RETURN
c
c  solve AX=B and AX=C
      DO 10 J=1,N
        CALL DGESL(A,M,M,WORKI,B(1,J),0)
   10 CONTINUE
      RETURN
      END






C***************************************************************
C
C     SUBROUTINE : MLINEAR2
C
C     SYNOPSIS :  Solves two linear systems, AX=B and
C                   AZ=C.  Condition number is not checked
C
C     PARAMETERS
C          A     ON INPUT:  MxM double precision matrix
C                ON OUTPUT: lu decompsoition of A
C          B     ON INPUT: MxN
C                ON OUTPUT  the solution X
C          C     ON INPUT:  MxN
C                ON OUTPUT  the solution Z
C
C***************************************************************

      SUBROUTINE MLINEAR2(A,B,C,M,N,WORKI,ERR)
      INTEGER*4 M,N, WORKI(M)
      DOUBLE PRECISION  A(M,M), B(M,N), C(M,N)
      LOGICAL ERR

c
c  local declarations
      INTEGER*4  INFO
c
c  get the lu decompsoition
      CALL DGEFA (A,M,M,WORKI,INFO)
      ERR = (INFO.NE.0)      
      IF (ERR) RETURN
c
c  solve AX=B and AX=C
      DO 10 J=1,N
        CALL DGESL(A,M,M,WORKI,B(1,J),0)
        CALL DGESL(A,M,M,WORKI,C(1,J),0)
   10 CONTINUE
      RETURN
      END


C***************************************************************
C
C     SUBROUTINE : VLINEAR
C
C     SYNOPSIS :  Solves the linear system, Ax=B when B is a vector
C
C     PARAMETERS
C          A     ON INPUT:  MxM double precision matrix
C                ON OUTPUT: lu decompsoition of A
C          B     ON INPUT: Mx1
C                ON OUTPUT  the solution x
C
C***************************************************************

      SUBROUTINE VLINEAR(A,B,M,WORKI,ERR)
      INTEGER*4 M, WORKI(M)
      DOUBLE PRECISION  A(M,M), B(M)
      LOGICAL ERR

c
c  local declarations
      INTEGER*4  INFO
c
c  get the lu decompsoition
      CALL DGEFA (A,M,M,WORKI,INFO)
      ERR = (INFO.NE.0)
      IF (ERR) RETURN
c
c  solve Ax=b
      CALL DGESL(A,M,M,WORKI,B,0)
c
c  finish
   10 CONTINUE
      RETURN
      END












C______________________________________________________________________
C
C   LINPACK  ROUTINES
C
C________________________________________________________________________


c THIS IS A MODIFIED LINPACK ROUTINE DESIGNED FOR 
c USE WITH MATLAB MEX FILES - THE ONLY CHANGE IS THAT INTEGER
c VARIABLES WERE CHANGED TO INTEGER*4 VARIABLES
c PROGRAM MODIFIED Nov 1994 by Evan Anderson
      subroutine dgedi(a,lda,n,ipvt,det,work,job)
      integer*4 lda,n,ipvt(1),job
      double precision a(lda,1),det(2),work(1)
c
c     dgedi computes the determinant and inverse of a matrix
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefa.
c
c        work    double precision(n)
c                work vector.  contents destroyed.
c
c        job     integer
c                = 11   both determinant and inverse.
c                = 01   inverse only.
c                = 10   determinant only.
c
c     on return
c
c        a       inverse of original matrix if requested.
c                otherwise unchanged.
c
c        det     double precision(2)
c                determinant of original matrix if requested.
c                otherwise not referenced.
c                determinant = det(1) * 10.0**det(2)
c                with  1.0 .le. dabs(det(1)) .lt. 10.0
c                or  det(1) .eq. 0.0 .
c
c     error condition
c
c        a division by zero will occur if the input factor contains
c        a zero on the diagonal and the inverse is requested.
c        it will not occur if the subroutines are called correctly
c        and if dgeco has set rcond .gt. 0.0 or dgefa has set
c        info .eq. 0 .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c THIS IS A MODIFIED LINPACK ROUTINE DESIGNED FOR 
c USE WITH MATLAB MEX FILES - THE ONLY CHANGE IS THAT INTEGER
c VARIABLES WERE CHANGED TO INTEGER*4 VARIABLES
c PROGRAM MODIFIED Nov 1994 by Evan Anderson
c
c     subroutines and functions
c
c     blas daxpy,dscal,dswap
c     fortran dabs,mod
c
c     internal variables
c
      double precision t
      double precision ten
      integer*4 i,j,k,kb,kp1,l,nm1
c
c
c     compute determinant
c
      if (job/10 .eq. 0) go to 70
         det(1) = 1.0d0
         det(2) = 0.0d0
         ten = 10.0d0
         do 50 i = 1, n
            if (ipvt(i) .ne. i) det(1) = -det(1)
            det(1) = a(i,i)*det(1)
c        ...exit
            if (det(1) .eq. 0.0d0) go to 60
   10       if (dabs(det(1)) .ge. 1.0d0) go to 20
               det(1) = ten*det(1)
               det(2) = det(2) - 1.0d0
            go to 10
   20       continue
   30       if (dabs(det(1)) .lt. ten) go to 40
               det(1) = det(1)/ten
               det(2) = det(2) + 1.0d0
            go to 30
   40       continue
   50    continue
   60    continue
   70 continue
c
c     compute inverse(u)
c
      if (mod(job,10) .eq. 0) go to 150
         do 100 k = 1, n
            a(k,k) = 1.0d0/a(k,k)
            t = -a(k,k)
            call dscal(k-1,t,a(1,k),1)
            kp1 = k + 1
            if (n .lt. kp1) go to 90
            do 80 j = kp1, n
               t = a(k,j)
               a(k,j) = 0.0d0
               call daxpy(k,t,a(1,k),1,a(1,j),1)
   80       continue
   90       continue
  100    continue
c
c        form inverse(u)*inverse(l)
c
         nm1 = n - 1
         if (nm1 .lt. 1) go to 140
         do 130 kb = 1, nm1
            k = n - kb
            kp1 = k + 1
            do 110 i = kp1, n
               work(i) = a(i,k)
               a(i,k) = 0.0d0
  110       continue
            do 120 j = kp1, n
               t = work(j)
               call daxpy(n,t,a(1,j),1,a(1,k),1)
  120       continue
            l = ipvt(k)
            if (l .ne. k) call dswap(n,a(1,k),1,a(1,l),1)
  130    continue
  140    continue
  150 continue
      return
      end





c THIS IS A MODIFIED LINPACK ROUTINE DESIGNED FOR 
c USE WITH MATLAB MEX FILES - THE ONLY CHANGE IS THAT INTEGER
c VARIABLES WERE CHANGED TO INTEGER*4 VARIABLES
c PROGRAM MODIFIED Nov 1994 by Evan Anderson
      subroutine dgefa(a,lda,n,ipvt,info)
      integer*4 lda,n,ipvt(1),info
      double precision a(lda,1)
c
c     dgefa factors a double precision matrix by gaussian elimination.
c
c     dgefa is usually called by dgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgesl or dgedi will divide by zero
c                     if called.  use  rcond  in dgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c THIS IS A MODIFIED LINPACK ROUTINE DESIGNED FOR 
c USE WITH MATLAB MEX FILES - THE ONLY CHANGE IS THAT INTEGER
c VARIABLES WERE CHANGED TO INTEGER*4 VARIABLES
c PROGRAM MODIFIED Nov 1994 by Evan Anderson
c
c     subroutines and functions
c
c     blas daxpy,dscal,idamax
c
c     internal variables
c
      double precision t
      integer*4 idamax,j,k,kp1,l,nm1
c
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0d0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0d0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end

c THIS IS A MODIFIED LINPACK ROUTINE DESIGNED FOR 
c USE WITH MATLAB MEX FILES - THE ONLY CHANGE IS THAT INTEGER
c VARIABLES WERE CHANGED TO INTEGER*4 VARIABLES
c PROGRAM MODIFIED Nov 1994 by Evan Anderson
      subroutine dgesl(a,lda,n,ipvt,b,job)
      integer*4 lda,n,ipvt(1),job
      double precision a(lda,1),b(1)
c
c     dgesl solves the double precision system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer*4
c                the leading dimension of the array  a .
c
c        n       integer*4
c                the order of the matrix  a .
c
c        ipvt    integer*4(n)
c                the pivot vector from dgeco or dgefa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c        job     integer*4
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if dgeco has set rcond .gt. 0.0
c        or dgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call dgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c THIS IS A MODIFIED LINPACK ROUTINE DESIGNED FOR 
c USE WITH MATLAB MEX FILES - THE ONLY CHANGE IS THAT INTEGER
c VARIABLES WERE CHANGED TO INTEGER*4 VARIABLES
c PROGRAM MODIFIED Nov 1994 by Evan Anderson
c
c     subroutines and functions
c
c     blas daxpy,ddot
c
c     internal variables
c
      double precision ddot,t
      integer*4 k,kb,l,nm1
c
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call daxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
c     
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call daxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
c     
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            t = ddot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
   60    continue
c     
c        now solve trans(l)*x = y
c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + ddot(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end











c THIS IS A MODIFIED LINPACK ROUTINE DESIGNED FOR 
c USE WITH MATLAB MEX FILES - THE ONLY CHANGE IS THAT INTEGER
c VARIABLES WERE CHANGED TO INTEGER*4 VARIABLES
c PROGRAM MODIFIED Nov 1994 by Evan Anderson

      subroutine dsidi(a,lda,n,kpvt,det,inert,work,job)
      integer*4 lda,n,job
      double precision a(lda,1),work(1)
      double precision det(2)
      integer*4 kpvt(1),inert(3)
c
c     dsidi computes the determinant, inertia and inverse
c     of a double precision symmetric matrix using the factors from
c     dsifa.
c
c     on entry
c
c        a       double precision(lda,n)
c                the output from dsifa.
c
c        lda     integer
c                the leading dimension of the array a.
c
c        n       integer
c                the order of the matrix a.
c
c        kpvt    integer(n)
c                the pivot vector from dsifa.
c
c        work    double precision(n)
c                work vector.  contents destroyed.
c
c        job     integer
c                job has the decimal expansion  abc  where
c                   if  c .ne. 0, the inverse is computed,
c                   if  b .ne. 0, the determinant is computed,
c                   if  a .ne. 0, the inertia is computed.
c
c                for example, job = 111  gives all three.
c
c     on return
c
c        variables not requested by job are not used.
c
c        a      contains the upper triangle of the inverse of
c               the original matrix.  the strict lower triangle
c               is never referenced.
c
c        det    double precision(2)
c               determinant of original matrix.
c               determinant = det(1) * 10.0**det(2)
c               with 1.0 .le. dabs(det(1)) .lt. 10.0
c               or det(1) = 0.0.
c
c        inert  integer(3)
c               the inertia of the original matrix.
c               inert(1)  =  number of positive eigenvalues.
c               inert(2)  =  number of negative eigenvalues.
c               inert(3)  =  number of zero eigenvalues.
c
c     error condition
c
c        a division by zero may occur if the inverse is requested
c        and  dsico  has set rcond .eq. 0.0
c        or  dsifa  has set  info .ne. 0 .
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab
c THIS IS A MODIFIED LINPACK ROUTINE DESIGNED FOR 
c USE WITH MATLAB MEX FILES - THE ONLY CHANGE IS THAT INTEGER
c VARIABLES WERE CHANGED TO INTEGER*4 VARIABLES
c PROGRAM MODIFIED Nov 1994 by Evan Anderson
c
c     subroutines and functions
c
c     blas daxpy,dcopy,ddot,dswap
c     fortran dabs,iabs,mod
c
c     internal variables.
c
      double precision akkp1,ddot,temp
      double precision ten,d,t,ak,akp1
      integer*4 j,jb,k,km1,ks,kstep
      logical noinv,nodet,noert
c
      noinv = mod(job,10) .eq. 0
      nodet = mod(job,100)/10 .eq. 0
      noert = mod(job,1000)/100 .eq. 0
c
      if (nodet .and. noert) go to 140
         if (noert) go to 10
            inert(1) = 0
            inert(2) = 0
            inert(3) = 0
   10    continue
         if (nodet) go to 20
            det(1) = 1.0d0
            det(2) = 0.0d0
            ten = 10.0d0
   20    continue
         t = 0.0d0
         do 130 k = 1, n
            d = a(k,k)
c
c           check if 1 by 1
c
            if (kpvt(k) .gt. 0) go to 50
c
c              2 by 2 block
c              use det (d  s)  =  (d/t * c - t) * t  ,  t = dabs(s)
c                      (s  c)
c              to avoid underflow/overflow troubles.
c              take two passes through scaling.  use  t  for flag.
c
               if (t .ne. 0.0d0) go to 30
                  t = dabs(a(k,k+1))
                  d = (d/t)*a(k+1,k+1) - t
               go to 40
   30          continue
                  d = t
                  t = 0.0d0
   40          continue
   50       continue
c
            if (noert) go to 60
               if (d .gt. 0.0d0) inert(1) = inert(1) + 1
               if (d .lt. 0.0d0) inert(2) = inert(2) + 1
               if (d .eq. 0.0d0) inert(3) = inert(3) + 1
   60       continue
c
            if (nodet) go to 120
               det(1) = d*det(1)
               if (det(1) .eq. 0.0d0) go to 110
   70             if (dabs(det(1)) .ge. 1.0d0) go to 80
                     det(1) = ten*det(1)
                     det(2) = det(2) - 1.0d0
                  go to 70
   80             continue
   90             if (dabs(det(1)) .lt. ten) go to 100
                     det(1) = det(1)/ten
                     det(2) = det(2) + 1.0d0
                  go to 90
  100             continue
  110          continue
  120       continue
  130    continue
  140 continue
c
c     compute inverse(a)
c
      if (noinv) go to 270
         k = 1
  150    if (k .gt. n) go to 260
            km1 = k - 1
            if (kpvt(k) .lt. 0) go to 180
c
c              1 by 1
c
               a(k,k) = 1.0d0/a(k,k)
               if (km1 .lt. 1) go to 170
                  call dcopy(km1,a(1,k),1,work,1)
                  do 160 j = 1, km1
                     a(j,k) = ddot(j,a(1,j),1,work,1)
                     call daxpy(j-1,work(j),a(1,j),1,a(1,k),1)
  160             continue
                  a(k,k) = a(k,k) + ddot(km1,work,1,a(1,k),1)
  170          continue
               kstep = 1
            go to 220
  180       continue
c
c              2 by 2
c
               t = dabs(a(k,k+1))
               ak = a(k,k)/t
               akp1 = a(k+1,k+1)/t
               akkp1 = a(k,k+1)/t
               d = t*(ak*akp1 - 1.0d0)
               a(k,k) = akp1/d
               a(k+1,k+1) = ak/d
               a(k,k+1) = -akkp1/d
               if (km1 .lt. 1) go to 210
                  call dcopy(km1,a(1,k+1),1,work,1)
                  do 190 j = 1, km1
                     a(j,k+1) = ddot(j,a(1,j),1,work,1)
                     call daxpy(j-1,work(j),a(1,j),1,a(1,k+1),1)
  190             continue
                  a(k+1,k+1) = a(k+1,k+1) + ddot(km1,work,1,a(1,k+1),1)
                  a(k,k+1) = a(k,k+1) + ddot(km1,a(1,k),1,a(1,k+1),1)
                  call dcopy(km1,a(1,k),1,work,1)
                  do 200 j = 1, km1
                     a(j,k) = ddot(j,a(1,j),1,work,1)
                     call daxpy(j-1,work(j),a(1,j),1,a(1,k),1)
  200             continue
                  a(k,k) = a(k,k) + ddot(km1,work,1,a(1,k),1)
  210          continue
               kstep = 2
  220       continue
c
c           swap
c
            ks = iabs(kpvt(k))
            if (ks .eq. k) go to 250
               call dswap(ks,a(1,ks),1,a(1,k),1)
               do 230 jb = ks, k
                  j = k + ks - jb
                  temp = a(j,k)
                  a(j,k) = a(ks,j)
                  a(ks,j) = temp
  230          continue
               if (kstep .eq. 1) go to 240
                  temp = a(ks,k+1)
                  a(ks,k+1) = a(k,k+1)
                  a(k,k+1) = temp
  240          continue
  250       continue
            k = k + kstep
         go to 150
  260    continue
  270 continue
      return
      end

c THIS IS A MODIFIED LINPACK ROUTINE DESIGNED FOR 
c USE WITH MATLAB MEX FILES - THE ONLY CHANGE IS THAT INTEGER
c VARIABLES WERE CHANGED TO INTEGER*4 VARIABLES
c PROGRAM MODIFIED Nov 1994 by Evan Anderson
      subroutine dsifa(a,lda,n,kpvt,info)
      integer*4 lda,n,kpvt(1),info
      double precision a(lda,1)
c
c     dsifa factors a double precision symmetric matrix by elimination
c     with symmetric pivoting.
c
c     to solve  a*x = b , follow dsifa by dsisl.
c     to compute  inverse(a)*c , follow dsifa by dsisl.
c     to compute  determinant(a) , follow dsifa by dsidi.
c     to compute  inertia(a) , follow dsifa by dsidi.
c     to compute  inverse(a) , follow dsifa by dsidi.
c
c     on entry
c
c        a       double precision(lda,n)
c                the symmetric matrix to be factored.
c                only the diagonal and upper triangle are used.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       a block diagonal matrix and the multipliers which
c                were used to obtain it.
c                the factorization can be written  a = u*d*trans(u)
c                where  u  is a product of permutation and unit
c                upper triangular matrices , trans(u) is the
c                transpose of  u , and  d  is block diagonal
c                with 1 by 1 and 2 by 2 blocks.
c
c        kpvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if the k-th pivot block is singular. this is
c                     not an error condition for this subroutine,
c                     but it does indicate that dsisl or dsidi may
c                     divide by zero if called.
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab.
c THIS IS A MODIFIED LINPACK ROUTINE DESIGNED FOR 
c USE WITH MATLAB MEX FILES - THE ONLY CHANGE IS THAT INTEGER
c VARIABLES WERE CHANGED TO INTEGER*4 VARIABLES
c PROGRAM MODIFIED Nov 1994 by Evan Anderson
c
c     subroutines and functions
c
c     blas daxpy,dswap,idamax
c     fortran dabs,dmax1,dsqrt
c
c     internal variables
c
      double precision ak,akm1,bk,bkm1,denom,mulk,mulkm1,t
      double precision absakk,alpha,colmax,rowmax
      integer*4 imax,imaxp1,j,jj,jmax,k,km1,km2,kstep,idamax
      logical swap
c
c
c     initialize
c
c     alpha is used in choosing pivot block size.
      alpha = (1.0d0 + dsqrt(17.0d0))/8.0d0
c
      info = 0
c
c     main loop on k, which goes from n to 1.
c
      k = n
   10 continue
c
c        leave the loop if k=0 or k=1.
c
c     ...exit
         if (k .eq. 0) go to 200
         if (k .gt. 1) go to 20
            kpvt(1) = 1
            if (a(1,1) .eq. 0.0d0) info = 1
c     ......exit
            go to 200
   20    continue
c
c        this section of code determines the kind of
c        elimination to be performed.  when it is completed,
c        kstep will be set to the size of the pivot block, and
c        swap will be set to .true. if an interchange is
c        required.
c
         km1 = k - 1
         absakk = dabs(a(k,k))
c
c        determine the largest off-diagonal element in
c        column k.
c
         imax = idamax(k-1,a(1,k),1)
         colmax = dabs(a(imax,k))
         if (absakk .lt. alpha*colmax) go to 30
            kstep = 1
            swap = .false.
         go to 90
   30    continue
c
c           determine the largest off-diagonal element in
c           row imax.
c
            rowmax = 0.0d0
            imaxp1 = imax + 1
            do 40 j = imaxp1, k
               rowmax = dmax1(rowmax,dabs(a(imax,j)))
   40       continue
            if (imax .eq. 1) go to 50
               jmax = idamax(imax-1,a(1,imax),1)
               rowmax = dmax1(rowmax,dabs(a(jmax,imax)))
   50       continue
            if (dabs(a(imax,imax)) .lt. alpha*rowmax) go to 60
               kstep = 1
               swap = .true.
            go to 80
   60       continue
            if (absakk .lt. alpha*colmax*(colmax/rowmax)) go to 70
               kstep = 1
               swap = .false.
            go to 80
   70       continue
               kstep = 2
               swap = imax .ne. km1
   80       continue
   90    continue
         if (dmax1(absakk,colmax) .ne. 0.0d0) go to 100
c
c           column k is zero.  set info and iterate the loop.
c
            kpvt(k) = k
            info = k
         go to 190
  100    continue
         if (kstep .eq. 2) go to 140
c
c           1 x 1 pivot block.
c
            if (.not.swap) go to 120
c
c              perform an interchange.
c
               call dswap(imax,a(1,imax),1,a(1,k),1)
               do 110 jj = imax, k
                  j = k + imax - jj
                  t = a(j,k)
                  a(j,k) = a(imax,j)
                  a(imax,j) = t
  110          continue
  120       continue
c
c           perform the elimination.
c
            do 130 jj = 1, km1
               j = k - jj
               mulk = -a(j,k)/a(k,k)
               t = mulk
               call daxpy(j,t,a(1,k),1,a(1,j),1)
               a(j,k) = mulk
  130       continue
c
c           set the pivot array.
c
            kpvt(k) = k
            if (swap) kpvt(k) = imax
         go to 190
  140    continue
c
c           2 x 2 pivot block.
c
            if (.not.swap) go to 160
c
c              perform an interchange.
c
               call dswap(imax,a(1,imax),1,a(1,k-1),1)
               do 150 jj = imax, km1
                  j = km1 + imax - jj
                  t = a(j,k-1)
                  a(j,k-1) = a(imax,j)
                  a(imax,j) = t
  150          continue
               t = a(k-1,k)
               a(k-1,k) = a(imax,k)
               a(imax,k) = t
  160       continue
c
c           perform the elimination.
c
            km2 = k - 2
            if (km2 .eq. 0) go to 180
               ak = a(k,k)/a(k-1,k)
               akm1 = a(k-1,k-1)/a(k-1,k)
               denom = 1.0d0 - ak*akm1
               do 170 jj = 1, km2
                  j = km1 - jj
                  bk = a(j,k)/a(k-1,k)
                  bkm1 = a(j,k-1)/a(k-1,k)
                  mulk = (akm1*bk - bkm1)/denom
                  mulkm1 = (ak*bkm1 - bk)/denom
                  t = mulk
                  call daxpy(j,t,a(1,k),1,a(1,j),1)
                  t = mulkm1
                  call daxpy(j,t,a(1,k-1),1,a(1,j),1)
                  a(j,k) = mulk
                  a(j,k-1) = mulkm1
  170          continue
  180       continue
c
c           set the pivot array.
c
            kpvt(k) = 1 - k
            if (swap) kpvt(k) = -imax
            kpvt(k-1) = kpvt(k)
  190    continue
         k = k - kstep
      go to 10
  200 continue
      return
      end

