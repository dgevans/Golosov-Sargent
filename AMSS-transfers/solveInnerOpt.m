function [n,xprime,c,exitflag] =solveInnerOpt(x,s_,coeff,V,Para)
g=Para.g;
psi=Para.psi;
beta=Para.beta;
xMax=Para.xMax;
xMin=Para.xMin;
options=Para.options;
pi=Para.pi;
sSize=Para.sSize;
nguess=[Para.nDet(x) Para.nDet(x)];
 get_root_labor_nag= @(num,n,user,iflag) getResLaborFsolve(num,n,x,s_,coeff,V,Para,user,iflag)  ;
 [n,~,exitflag]=c05qb(get_root_labor_nag,nguess,'xtol',1e-10);
%[n,~,exitflag]=fsolve(@(n)getResLaborFsolve(x,s_,n,coeff,V,Para),nguess ,options);
c=n-g;
R=1./(beta*(c).*sum(pi(s_,:).*(1./c)));
% Use Implementability to get xprime
xprime=(n./(1-n))*(1-psi)+ x.*R-psi;

%% Check bounds
flagConsBind=0;
for s=1:sSize
if xprime(s)>xMax
    xprime(s)=xMax;
    flagConsBind=1;
elseif  xprime(s)<xMin
     xprime(s)=xMin;
    flagConsBind=1;
end
end

if flagConsBind>0
[n,~,exitflag,message,jacob]=fsolve(@(n)getResInitialLaborFsolve(x,s_,n,Para,squeeze(xprime)'),n,options);
c=n-g;
R=1./(beta*(c).*sum(pi(s_,:).*(1./c)));
% Use Implementability to get xprime
xprime=(n./(1-n))*(1-psi)+ x.*R-psi;
if exitflag==1
    exitflag=0;
end
end


end
