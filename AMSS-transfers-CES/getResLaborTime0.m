function [resEQ0,user,iflag]=getResLaborTime0(num,b_,user,n0,coeff,V,Para,iflag)
pi=Para.pi;
g=Para.g;
beta=Para.beta;
sSize=Para.sSize;
der_u_c=Para.der_u_c;
der_u_cn=Para.der_u_cn;
der_u_n=Para.der_u_n;
der_u_nn=Para.der_u_nn;

eq=[];
eqGrad=[];
u_c=der_u_c(n0);
u_n=der_u_n(n0);
u_nn=der_u_nn(n0);
u_cn=der_u_cn(n0);
Euc=sum(pi(s_,:).*u_c);
R=u_c./(beta*Euc);

xprime0=u_n.*n-(n-g).*u_c + u_c*b_;

for s=1:sSize
VJac(s)=(u_c(s)-u_n(s))*pi(s_,s);
phi(s)=beta*funeval(coeff(s,:)',V(s),xprime0(s),1)*pi(s_,s);
end

% Building the Jacobian of I
WithPostTaxLaborIncome =diag(u_n+n.*u_nn);
With_b_terms=(b_)*(diag(u_cn));

ImpConsJac=WithPostTaxLaborIncome+With_b_terms;

resEQ=VJac+phi*ImpConsJac;

