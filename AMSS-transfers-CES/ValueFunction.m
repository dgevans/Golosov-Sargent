%function [mode, VNew, VJac, user] = ValueFunction(mode, n, z,s_,coeff,V,Para,objgrd, nstate, user)
function [VNew, VJac] = ValueFunction( z,s_,coeff,V,Para)

pi=Para.pi;
beta=Para.beta;
sSize=Para.sSize;
der_u_n=Para.der_u_n;
util=Para.util;
n=z(1:sSize);
xprime=z(sSize+1:end);

u_n=Para.der_u_n(n);
u_c=Para.der_u_c(n);
for s=1:sSize
VJac(s)=(u_c(s)-u_n(s))*pi(s_,s);
VJac(sSize+s)=beta*funeval(coeff(s,:)',V(s),xprime(s),1)*pi(s_,s);
Vstar(s)=funeval(coeff(s,:)',V(s),xprime(s));
end
u=util(n)+beta*Vstar;
VNew=-sum(pi(s_,:).*u);
VJac=-VJac;
end
