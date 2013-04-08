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
end
u=psi*log((n-g))+(1-psi)*log(1-n)+beta*Vstar;
VNew=-sum(pi(s_,:).*u);
VJac(1)=(psi/(n(1)-g(1))-(1-psi)/(1-n(1)))*pi(s_,1);
VJac(2)=(psi/(n(2)-g(2))-(1-psi)/(1-n(2)))*pi(s_,2);
VJac(3)=beta*funeval(coeff(1,:)',V(1),xprime(1),1)*pi(s_,1);
VJac(4)=beta*funeval(coeff(2,:)',V(2),xprime(2),1)*pi(s_,2);
VJac=-VJac;
end
