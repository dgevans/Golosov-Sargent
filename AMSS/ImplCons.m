function [ResINQ,ResEQ]=ImplCons(yy,x,s_,coeff,V,Para)

g=Para.g;
psi=Para.psi;
pi=Para.pi;
beta=Para.beta;
sSize=Para.sSize;
n=yy(1:sSize);
xprime=yy(sSize+1:end);
c=n-g;
R=1./(beta*(c).*sum(pi(s_,:).*(1./c)));
ResEQ=(n./(1-n))*(1-psi)+ x.*R-psi-xprime;

ResINQ=[];