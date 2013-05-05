%function [mode, ImpCons, ImpConsJac, user] = ImplementabilityCons(mode, ncnln, n, ldcj, needc, z, x,s_,Para,ImpConsJac, nstate, user)
function [ImpCons,eq,ImpConsJac,eqGrad] = ImplementabilityCons(z, x,s_,Para)
g=Para.g;
psi=Para.psi;
pi=Para.pi;
beta=Para.beta;
sSize=Para.sSize;
eq=[];
eqGrad=[];
n=z(1:sSize);
xprime=z(sSize+1:end);
c=n-g;
R=1./(beta*(c).*sum(pi(s_,:).*(1./c)));  
u_c=1./c;
u_cn=-1./(c.^2);
Euc=sum(pi(s_,:).*u_c);
ImpCons=x.*R+(1-psi)*(n./(1-n))-psi-xprime;
S=Para.sSize;

% Building the Jacobian of I
WithPostTaxLaborIncome =diag((1-psi)./(1-n).^2);
With_x_terms=(x/beta)*(diag(u_cn./Euc)-kron(u_c',u_cn.*pi(s_,:)./Euc^2));

With_xprime =-eye(S);

ImpConsJac=horzcat(WithPostTaxLaborIncome+With_x_terms ,With_xprime)';
