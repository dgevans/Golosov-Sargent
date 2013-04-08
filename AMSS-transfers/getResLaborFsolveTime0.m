function [resEQ0]=getResLaborFsolveTime0(b_,s0,n0,coeff,V,Para)

g=Para.g;
psi=Para.psi;
pi=Para.pi;
beta=Para.beta;
sSize=Para.sSize;
c0=n0-g(s0);
% Implementability at s0 and b_
xprime0=(n0/(1-n0))*(1-psi)+ (psi/c0)*b_-psi;
if min(n0-g(s0))>0 && min(1-n(0))>0
phi0=beta*funeval(coeff(s0,:)',V(s0),xprime0,1);
term1=(psi/(n0-g(s0))-(1-psi)/(1-n0));
term2=-(1-psi)/(1-n0)^2;
resEQ0=term1-phi0*term2;
else
  resEQ0=(abs(n0-g(s0))+abs(1-n0)).*10;
end

