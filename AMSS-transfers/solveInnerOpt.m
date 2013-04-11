function [n,xprime,c,VNew,exitflag] =solveInnerOpt(x,s_,coeff,V,Para)
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
 [n,~,exitflag]=c05qb(get_root_labor_nag,nguess,'xtol',1e-13);
 
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

% if Para.flagTransfers==2
% for s=1:sSize    
% phi(s)=beta*funeval(coeff(s,:)',V(s),xprime(s),1);
% if phi(s)>-1e-4
%     n=Para.n_fb;
%     
% c=n-g;
% % Use Implementability to get xprime
% xprime=Para.x_fb(1)*ones(1,sSize);
%     
% end

%end
%end
for s=1:sSize
Vstar(s)=funeval(coeff(s,:)',V(s),xprime(s));
end
u=psi*log((n-g))+(1-psi)*log(1-n)+beta*Vstar;
VNew=sum(pi(s_,:).*u);

end
