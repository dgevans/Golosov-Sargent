
function [res]=resCrossingPoints(x,coeff,V,Para)
[~,xprime,~,~] =solveInnerOpt(x,1,coeff,V,Para);
res=xprime(1)-xprime(2);
end