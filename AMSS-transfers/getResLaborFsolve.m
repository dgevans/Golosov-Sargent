function [resEQ,user,iflag]=getResLaborFsolve(num,n,x,s_,coeff,V,Para,user,iflag)

g=Para.g;
psi=Para.psi;
pi=Para.pi;
beta=Para.beta;
sSize=Para.sSize;
c=n-g;
R=1./(beta*(c).*sum(pi(s_,:).*(1./c)));  
u_c=1./c;
u_cn=-1./(c.^2);
Euc=sum(pi(s_,:).*u_c);
S=Para.sSize;

if min(n-g)>0 && min(1-n)>0
    % Use Implementability to get xprime
xprime=(n./(1-n))*(1-psi)+ x.*R-psi;

for s=1:sSize
VJac(s)=(psi/(n(s)-g(s))-(1-psi)/(1-n(s)))*pi(s_,s);
phi(s)=beta*funeval(coeff(s,:)',V(s),xprime(s),1)*pi(s_,s);
end

% Building the Jacobian of I
WithPostTaxLaborIncome =diag((1-psi)./(1-n).^2);
With_x_terms=(x/beta)*(diag(u_cn./Euc)-kron(u_c',u_cn.*pi(s_,:)./Euc^2));

ImpConsJac=WithPostTaxLaborIncome+With_x_terms;

resEQ=VJac+phi*ImpConsJac;

else
  resEQ=(abs(n-g)+abs(1-n)).*10;
end


