% Inputs - xInit, state variables - u2btild,,R,s_  coeff, value
% function, para
function [V_new V_other]=CheckOpt(u2bdiff,RR,s,c,VV,xInit,Para)
global V Vcoef R u2btild Par s_ flagCons

%Get the initial guess for the uconstraint problem. With the simplification
%we need only c1_1,c1_2and c2_1

xInit=xInit(1:3);

errValue=0;
Para.theta=[Para.theta_1 Para.theta_2];
Para.alpha=[Para.alpha_1 Para.alpha_2];
Par=Para;
u2btild=u2bdiff;
R=RR;
Vcoef{1}=c(1,:)';
Vcoef{2}=c(2,:)';
V=VV;
s_=s;


f = @(x) -BelObjectiveUncond(3,x,1);


[~,V_new] = fminsearch(f,xInit);
V_new = -V_new;
V_other = -f(xInit);
end
