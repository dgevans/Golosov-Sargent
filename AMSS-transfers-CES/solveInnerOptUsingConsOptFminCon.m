function [n,xprime,c,VNew,exitflag] =solveInnerOptUsingConsOptFminCon(x,s_,coeff,V,Para,z)
g=Para.g;
sSize=Para.sSize;
a = [];
bl = vertcat(Para.g',Para.x_fb*ones(sSize,1));
bu = vertcat(Inf*ones(sSize,1),Para.xMax*ones(sSize,1));
confun=@(z) ImplementabilityCons(z, x,s_,Para);
objfun=@(z) ValueFunction(z,s_,coeff,V,Para);
options = optimset('Algorithm','interior-point');
options = optimset(options,'GradObj','on','GradConstr','on','Display','off'); 
%[z,VNew,exitflag,~,lambda] = fmincon(objfun,z,[],[],[],[],bl,bu,confun,options);

[z,VNew,exitflag,~,lambda] = ktrlink(objfun,z,[],[],[],[],bl,bu,confun,options);
 n=z(1:sSize);
xprime=z(sSize+1:end);
if exitflag>0
    exitflag=Para.solveflag;
end
VNew=-VNew;
c=n-g;
end
