function [ gapprox ] = fgradapprox( f,x,h )
%FGRADAPPROX Summary of this function goes here
%   Detailed explanation goes here
F = f(x);
n = length(F);
m = length(x);
xtild = x;
gapprox = zeros(m,n);
for i = 1:m
    xtild(i) = x(i) +h;
    Ftild = f(xtild);
    for j = 1:n
        gapprox(i,j) = (Ftild(j)-F(j))/h;
    end
    xtild = x;
end

end

