function [resEQ]=getResInitialLaborFsolve(x,s_,n,Para,xprime)
pi=Para.pi;
beta=Para.beta;
g=Para.g;
der_u_c=Para.der_u_c;
der_u_n=Para.der_u_n;
u_c=der_u_c(n);
u_n=der_u_n(n);
Euc=sum(pi(s_,:).*u_c);
R=u_c./(beta*Euc);

    % Use Implementability to get xprime
if min(n-g)>0 
resEQ=-xprime+u_n.*n-(n-g).*u_c + x.*R;

else
    resEQ=(abs(n-g)+abs(1-n)).*10;
end
