function [Tau0,Rprime0,u2btildprime0]=SolveTime0(c,V,s_,Para)
n1=Para.n1;
n2=Para.n2;
alpha_1=Para.alpha_1;
alpha_2=Para.alpha_2;
g=Para.g(s_);
theta_1=Para.theta_1;
theta_2=Para.theta_2;
psi=Para.psi;
beta=Para.beta;
% simulation
btild_1=Para.btild_1;
  disp('Computed V...Now solving V0(btild_1) where btild_1 is')
disp(btild_1)
% c1 and c2 solve 
 options=optimset('Display','off');
[x,fval,exitflagv0,~,grad] = fminunc(@(x)  getValue0(x, btild_1,1,Para,c,V),[ 1 mean(Para.RGrid)^(-1)],options);
if ~(exitflagv0==1)
    disp('Optimization failed for V0 once ..trying with fmincon')
opts = optimset('Algorithm', 'interior-point', 'Display','off', ...
    'GradObj','off','GradConstr','off',...
    'MaxIter',1000, ...
    'TolX', Para.ctol/10, 'TolFun', Para.ctol, 'TolCon', Para.ctol,'MaxTime',200);
lb=[0.001 0.001];
ub=[10 10];
%[x,fval,exitflagv0,output,lambda]  =fmincon(@(x) getValue0(x, btild_1,1,Para,c,V),[ x ],[],[],[],[],lb,ub,[],opts);
[x,fval,exitflagv0,output,lambda]  =fmincon(@(x) getValue0(x, btild_1,1,Para,c,V),[ 1 mean(Para.RGrid)^(-1)],[],[],[],[],lb,ub,[],opts);
end
c10 = x(1);
c20 = x(2);
R0=c10/c20;
TotalResources=(c10*n1+c20*n2+g(1));
FF=R0*theta_2/theta_1;
DenL2=n1*theta_1*FF+theta_2*n2;
l20=(TotalResources-n1*theta_1+n1*theta_1*FF)/(DenL2);
l10= 1-FF*(1-l20);
BracketTerm=l20/(1-l20)-(l10/(1-l10))*R0;
u2btildprime0=(((1-psi)/(psi))*BracketTerm+btild_1/(beta*psi)+R0-1)*psi;
btildprime0=u2btildprime0/(c20^-1*psi) ;
Rprime0=c20^(-1)/c10^(-1);
ul0=(1-psi)/(1-l10);
uc0=psi/c10;

Tau0=1-(ul0/(theta_1*uc0));

end
