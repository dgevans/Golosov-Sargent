function [n,xprime,c,VNew,exitflag] =solveInnerOptUsingConsOptFminCon(x,s_,coeff,V,Para,z)
g=Para.g;
psi=Para.psi;
beta=Para.beta;
xMax=Para.xMax;
xMin=Para.xMin;
options=Para.options;
pi=Para.pi;
sSize=Para.sSize;
nguess=[Para.nDet(x)]*ones(1,sSize);
c=nguess-g;
R=1./(beta*(c).*sum(pi(s_,:).*(1./c)));
% Use Implementability to get xprime
a = [];
bl = vertcat(Para.g',Para.x_fb*ones(sSize,1));
bu = vertcat(ones(sSize,1),Para.xMax*ones(sSize,1));
confun=@(z) ImplementabilityCons(z, x,s_,Para);
objfun=@(z) ValueFunction(z,s_,coeff,V,Para);
options = optimset('Algorithm','active-set');
options = optimset(options,'GradObj','on','GradConstr','on','Display','off');
[z,VNew,exitflag,~,lambda] = fmincon(objfun,z,[],[],[],[],bl,bu,... 
   confun,options);
n=z(1:sSize);
xprime=z(sSize+1:end);
c=n-g;
if exitflag>0
    exitflag=Para.solveflag;
end
VNew=-VNew;
end
