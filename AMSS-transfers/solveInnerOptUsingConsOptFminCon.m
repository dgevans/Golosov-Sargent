function [n,xprime,c,VNew,exitflag] =solveInnerOptUsingConsOptFminCon(x,s_,coeff,V,Para,z)
g=Para.g;
psi=Para.psi;
beta=Para.beta;
xMax=Para.xMax;
xMin=Para.xMin;
options=Para.options;
pi=Para.pi;
sSize=Para.sSize;
nguess=[Para.nDet(x) Para.nDet(x)];
c=nguess-g;
R=1./(beta*(c).*sum(pi(s_,:).*(1./c)));
% Use Implementability to get xprime
xprimeguess=(nguess./(1-nguess))*(1-psi)+ x.*R-psi;
a = [];
bl = [Para.g(1); Para.g(2); Para.xMin; Para.xMin];
bu = [1; 1; Para.xMax; Para.xMax];
confun=@(z) ImplementabilityCons(z, x,s_,Para);
objfun=@(z) ValueFunction(z,s_,coeff,V,Para);
options = optimset('Algorithm','interior-point');
options = optimset(options,'GradObj','on','GradConstr','on','Display','off');
[z,VNew,exitflag,~,lambda] = fmincon(objfun,z,[],[],[],[],bl,bu,... 
   confun,options);
n=z(1:2);
xprime=z(3:4);
c=n-g;
if exitflag>0
    exitflag=1;
end
VNew=-VNew;
end
