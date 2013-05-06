function [n,xprime,c,VNew,exitflag] =solveInnerOpt(x,s_,coeff,V,Para,nguess)
g=Para.g;
der_u_c=Para.der_u_c;
der_u_n=Para.der_u_n;
util=Para.util;
beta=Para.beta;
xMax=Para.xMax;
xMin=Para.xMin;
options=Para.options;
pi=Para.pi;
sSize=Para.sSize;
 get_root_labor_nag= @(num,n,user,iflag) getResLaborFsolve(num,n,x,s_,coeff,V,Para,user,iflag)  ;
 [n,~,exitflag]=c05qb(get_root_labor_nag,nguess);
 if ~(exitflag==0)
     disp(x)
 end
u_c=der_u_c(n); % marginal utility of consumption
Euc=sum(pi(s_,:).*u_c); % expected marginal utility of consumption
R=u_c./(beta*Euc); % factor that multiplies x
u_n=der_u_n(n); % marginal disutility of labor

% Use Implementability to get xprime
xprime=u_n.*n-(n-g).*u_c + x.*R;
flagConsBind=1;
%% Check bounds
n_check=1;
while (flagConsBind==1 )&& (n_check < sSize)

% check xprime    
for s=1:sSize
if xprime(s)>xMax
    xprime(s)=xMax;
    flagConsBind=1;
elseif  xprime(s)<xMin
     xprime(s)=xMin;
    flagConsBind=1;
else
    flagConsBind=0;
    
end
end

if flagConsBind>0
[n,~,exitflag,~,~]=fsolve(@(n)getResInitialLaborFsolve(x,s_,n,Para,xprime),n,options);


u_c=der_u_c(n);

end
n_check=n_check+1;
end

if exitflag==1
    exitflag=Para.solveflag;
end

Vstar=ones(1,sSize);

for s=1:sSize
Vstar(s)=funeval(coeff(s,:)',V(s),xprime(s));
end
u=util(n)+beta*Vstar;
VNew=sum(pi(s_,:).*u);
c=n-g;
