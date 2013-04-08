
function nag_roots_sys_func_easy_example 
% The following starting values provide a rough solution.
x = -ones(9, 1);
dd=2;
ff=@(n,x,user,iflag) fcn(n,x,user,dd,iflag)
[xOut, fvec, user, ifail] = nag_roots_sys_func_easy(ff, x);
switch ifail
  case {0}
    fprintf('\nFinal 2-norm of the residuals = %12.4e\n', norm(fvec));
    fprintf('\nFinal approximate solution\n');
    disp(xOut);
  case {2, 3, 4}
    fprintf('\nApproximate solution\n');
    disp(xOut);
end

end


function [fvec, user, iflag] = fcn(n, x, user, dd,iflag)
  nd = double(n); % n is an int64 and can't be used everywhere
  fvec = zeros(nd, 1);
  fvec(1:nd) = (3.0-2.0.*x).*x + 1.0*dd;
  fvec(2:nd) = fvec(2:nd) - x(1:(nd-1));
  fvec(1:(nd-1)) = fvec(1:(nd-1)) - 2.0.*x(2:nd);

end
 

