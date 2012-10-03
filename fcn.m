function [fvec, iflag] = fcn05nb(n,x,iflag)
  fvec = zeros(n, 1);
  for k = 1:double(n) 
    fvec(k) = (3.0-2.0*x(k))*x(k)+1.0;
    if k > 1
      fvec(k) = fvec(k) - x(k-1);
    end
    if k < n
      fvec(k) = fvec(k) - 2*x(k+1);
    end
  end