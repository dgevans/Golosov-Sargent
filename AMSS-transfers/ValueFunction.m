%function [mode, VNew, VJac, user] = ValueFunction(mode, n, z,s_,coeff,V,Para,objgrd, nstate, user)
function [VNew, VJac] = ValueFunction( z,s_,coeff,V,Para)

g=Para.g;
psi=Para.psi;
pi=Para.pi;
beta=Para.beta;
sSize=Para.sSize;

n=z(1:sSize);
xprime=z(sSize+1:end);

for s=1:sSize
Vstar(s)=funeval(coeff(s,:)',V(s),xprime(s));
VJac(s)=(psi/(n(s)-g(s))-(1-psi)/(1-n(s)))*pi(s_,s);
VJac(sSize+s)=beta*funeval(coeff(s,:)',V(s),xprime(s),1)*pi(s_,s);
end
u=psi*log((n-g))+(1-psi)*log(1-n)+beta*Vstar;
VNew=-sum(pi(s_,:).*u);
VJac=-VJac;
end
