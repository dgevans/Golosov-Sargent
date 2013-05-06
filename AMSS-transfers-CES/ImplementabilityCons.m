%function [mode, ImpCons, ImpConsJac, user] = ImplementabilityCons(mode, ncnln, n, ldcj, needc, z, x,s_,Para,ImpConsJac, nstate, user)
function [ImpCons,eq,ImpConsJac,eqGrad] = ImplementabilityCons(z, x,s_,Para)
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
n=z(1:sSize);
xprime=z(sSize+1:end);
u_c=der_u_c(n);
u_n=der_u_n(n);
u_nn=der_u_nn(n);
u_cn=der_u_cn(n);
Euc=sum(pi(s_,:).*u_c);
R=u_c./(beta*Euc);

ImpCons=u_n.*n-(n-g).*u_c + x.*R-xprime;
S=Para.sSize;

% Building the Jacobian of I
WithPostTaxLaborIncome =diag(u_n+n.*u_nn);
With_x_terms=(x/beta)*(diag(u_cn./Euc)-kron(u_c',u_cn.*pi(s_,:)./Euc^2));
With_xprime =-eye(S);
ImpConsJac=horzcat(WithPostTaxLaborIncome+With_x_terms ,With_xprime)';
