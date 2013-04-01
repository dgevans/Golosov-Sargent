function [V, LBAR, NBAR,alpha,beta,info,W]  = schurg(L,N,order,algorithm,eps)
%
% SCHURG  Ordered real generalized Schur decomposition
%
%
%  [V, LBAR,NBAR,alpha,beta,info,W]] = schurg(L,N) computes the ordered generalized real
%  Schur decomposition of the matrix pencil  
%                     lambda L  - N
%  such that LBAR is upper triangular, NBAR is upper block
%  triangular, V is the matrix of right Schur vectors such that
%  for some orthogonal matrix W
%
%          W' L V = LBAR
%          W' N V = NBAR
%
%   and the generalized eigenvalues of the pencil are given
%   by  (alpha ./ beta). 
%
%  Order, algorithm, and eps are optional parameters.  If they are
%  not present then they are given default values. 
%
%       If order =0 then the stable eigenvalues appear first in the pencil
%       lambda LBAR - NBAR. If order =1 then the unstable eigenvalues  appear first.
%       The default value of order is zero.
%
%       If algorithm = 0 then the decomposition is computed using
%       QZHESW, QZITW, QZVAL, and ORDER from  RICPACK. In this
%       algorithm the problem is not balanced.  If algorithm =1
%       then the same algorithm is used except that the problem is
%       balanced using   BALGEN and  BALGBK.  If algorithm
%       =2 then the decomposition is computed using DGGES from LAPACK.
%       The default value of algorithm is two.  
%
%       eps is a real number reflecting when an element of a matrix should
%       be considered zero. 
%
%
%  NOTES:
%
%    (1) If algorithm is equal to zero or one then the order of the
%    eigenvalues in 
%
%              (alpha ./ beta)
%
%     are not related to the order of the eigenvalues in the pencil.
%
%    (2) eps is only used when algorithm zero or one is selected.
%
%    (3) When algorithm two is used the matrix W is computed.
%    When algorithm zero or one is used (for efficiency reasons) it is not
%    computed.  In this case W is returned as a matrix of zeros.
%
%    (4)  L and N must be real matrices.  If they are
%             complex, the complex portion will be
%             truncated without warning.
%
%    (5) The output info is a two-dimensional vector giving information
%     on failure when either algorithm zero or one is used.  The casual user
%     need not worry about this since the algorithm almost never
%     fails.  If info(1) ~=0 then the algorithm failed to compute
%    the ordered generalized real Schur decomposition. If info(1)=1 this
%    is because the FORTRAN function QZITW failed. If info(1) =2 
%    this is because the FORTRAN function ORDER failed. Currently,
%    info(2) does not contain useful information.
%
%
%
% EXAMPLE:
%
%   (1)  To obtain a generalized real Schur decomposition
%     of the pencil lambda L  - N in which the stable
%     eigenvalues appear first, use the following code:
%            [V,LBAR,NBAR] = schurg(L,N)
%     If half of the eigenvalues of the pencil are stable and
%     if the algorithm didn't fail, the stabilizing solution to the
%     corresponding  Riccati like equation is
%          P = V21/V11
%     where V21 is the lower left block of V and V11 is the upper left
%     block of V.




%
%  HISTORY:
%
%  Versions:
%
%    1.0 :   December, 1994    Initial version 
%    1.1 :   January,1995      First general release. Documentation
%                                 was improved. 
%    1.2 :   January, 1999     First version to support windows
%    1.3 :   June, 2000        First version to support 64-bit
%                                    machines (alpha, SGI)
%    2.0 :   April, 2001       Added support for LAPACK
%    2.1 :   March 13, 2002      Added W as a possible output
%    3.0 :   November 4, 2002    Interface was made more user friendly  
%    3.0a:   June 19, 2003       Documentation improved.  
%                           






% Create the inputs,  if necessary.

if nargin< 3,
  order =0;
end
if nargin < 4,
  algorithm =2;
end
if nargin < 5,
  eps = 1e-15;
end


% Call the computational routine 

[V, LBAR, NBAR,alpha,beta,info,W]  = schurgaux(L,N,order,algorithm,eps);


