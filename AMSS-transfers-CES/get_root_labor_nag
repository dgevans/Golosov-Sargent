function [resEQ,iflag]=getResLaborFsolve(n,x,s_,coeff,V,Para,iflag)

g=Para.g;
psi=Para.psi;
pi=Para.pi;
beta=Para.beta;
sSize=Para.sSize;
c=n-g;
R=1./(beta*(c).*sum(pi(s_,:).*(1./c)));

% Use Implementability to get xprime
xprime=(n./(1-n))*(1-psi)+ x.*R-psi;
if min(n-g)>0 && min(1-n)>0

for s=1:sSize
    not_s=2/s;
phi(s)=beta*funeval(coeff(s,:)',V(s),xprime(s),1)*pi(s_,s);
phi(not_s)=beta*funeval(coeff(not_s,:)',V(not_s),xprime(not_s),1)*pi(s_,not_s);
term1(s)=(psi/(n(s)-g(s))-(1-psi)/(1-n(s)))*pi(s_,s);
term2a(s)=(1-psi)/(1-n(s))^2;
term2b(s)=x/(beta*(c(s))^2*sum(pi(s_,:).*(1./c)));
term2c(s)=pi(s_,s)/(c(s)*(sum(pi(s_,:).*(1./c))));
term3a(s)=phi(not_s);
term3b(s)=x/(beta*(n(not_s)-g(not_s)));
term3c(s)=(sum(pi(s_,:).*1./c))^(-2);
term3d(s)=pi(s_,s)/((n(s)-g(s))^2);
resEQ(s)=term1(s)+phi(s)*(term2a(s)+term2b(s)*(term2c(s)-1))+term3a(s)*term3b(s)*term3c(s)*term3d(s);
%resEQ(s)=term1(s)/((phi(s)*phi(not_s)))+(term2a(s)+term2b(s)*(term2c(s)-1))/(phi(not_s))+term3b(s)*term3c(s)*term3d(s)/(phi(s));
end
else
  resEQ=(abs(n-g)+abs(1-n)).*10;
end
