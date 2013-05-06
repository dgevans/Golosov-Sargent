function [DerL]=checkGradients(x,s_,n,coeff,V,Para,delta)
g=Para.g;
psi=Para.psi;
pi=Para.pi;
beta=Para.beta;
sSize=Para.sSize;
c=n-g;
R=1./(beta*(c).*sum(pi(s_,:).*(1./c)));

xprime=(n./(1-n))*(1-psi)+ x.*R-psi;
for s=1:2
phi(s)=beta*funeval(coeff(s,:)',V(s),xprime(s),1)*pi(s_,s);
end
pert=[delta 0];
DerL(1)=(Lagrangian(x,s_,n+pert,coeff,V,Para,phi,1)-Lagrangian(x,s_,n,coeff,V,Para,phi,1))/(delta);
pert=[0 delta];
DerL(2)=(Lagrangian(x,s_,n+pert,coeff,V,Para,phi,2)-Lagrangian(x,s_,n,coeff,V,Para,phi,2))/(delta);
DerL=DerL;
end


function L=Lagrangian(x,s_,n,coeff,V,Para,phi,ss)
g=Para.g;
psi=Para.psi;
pi=Para.pi;
beta=Para.beta;
sSize=Para.sSize;
c=n-g;
R=1./(beta*(c).*sum(pi(s_,:).*(1./c)));
L=0;
for s=1:2
Imp(s)=(-(n(s)./(1-n(s)))*(1-psi)- x.*R(s));    
LL(s)=(psi*log(n(s)-g(s))+(1-psi)*log(1-n(s)))*pi(s_,s)-phi(s)*Imp(s);
end
L=sum(LL);
end





