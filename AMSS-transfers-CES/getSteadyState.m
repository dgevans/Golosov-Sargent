function [res,user,iflag]=getSteadyState(num,ssPol,Para,user,iflag)
s_=1;
der_u_c=Para.der_u_c;
der_u_cn=Para.der_u_cn;
der_u_n=Para.der_u_n;
der_u_nn=Para.der_u_nn;
g=Para.g;
pi=Para.pi;
beta=Para.beta;
sSize=Para.sSize;
n(1)=ssPol(1);
n(2)=ssPol(2);
phiss=ssPol(3);
x=ssPol(4);

u_c=der_u_c(n);
u_n=der_u_n(n);
u_nn=der_u_nn(n);
u_cn=der_u_cn(n);
Euc=sum(pi(s_,:).*u_c);
R=u_c./(beta*Euc);


if min(n-g)>0 
    % Use Implementability to get xprime
xprime=u_n.*n-(n-g).*u_c + x.*R;

for s=1:sSize
VJac(s)=(u_c(s)-u_n(s))*pi(s_,s);
phi(s)=phiss;
end

% Building the Jacobian of I
WithPostTaxLaborIncome =diag(u_n+n.*u_nn);
With_x_terms=(x/beta)*(diag(u_cn./Euc)-kron(u_c',u_cn.*pi(s_,:)./Euc^2));

ImpConsJac=WithPostTaxLaborIncome+With_x_terms;

resEQ=VJac+phi*ImpConsJac;

else
  resEQ=(abs(n-g)).*10;
end


res=[resEQ xprime(1)-x xprime(2)-x];
