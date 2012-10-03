function [fvec, iflag] = fcn05nb(n,x,iflag)
   nd = double(n); % n is an int64 and can't be used everywhere
  fvec = zeros(nd, 1);
  fvec(1:nd) = (3.0-2.0.*x).*x + 1.0;
  fvec(2:nd) = fvec(2:nd) - x(1:(nd-1));
  fvec(1:(nd-1)) = fvec(1:(nd-1)) - 2.0.*x(2:nd);
  end