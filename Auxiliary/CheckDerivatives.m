% Inputs - xInit, state variables - u2btild,,R,s_  coeff, value
% function, para
function [DiffDerivatives]=CheckDerivatives(u2bdiff,RR,s,c,VV,x,Para,h)
global V Vcoef R u2btild Par s_

%Get the initial guess for the uconstraint problem. With the simplification
%we need only c1_1,c1_2and c2_1
Para.theta=[Para.theta_1 Para.theta_2];
Para.alpha=[Para.alpha_1 Para.alpha_2];
Par=Para;
u2btild=u2bdiff;
R=RR;
Vcoef{1}=c(1,:)';
Vcoef{2}=c(2,:)';
V=VV;
s_=s;
u2btildLL=Para.u2btildLL;
u2btildUL=Para.u2btildUL;
n1=Para.n1;
n2=Para.n2;
ctol=Para.ctol;

Der(1)= (Value3cont([x(1)+h,x(2),x(3)])-Value3cont([x(1)-h,x(2),x(3)]))/(2*h);
Der(2)= (Value3cont([x(1),x(2)+h,x(3)])-Value3cont([x(1),x(2)-h,x(3)]))/(2*h);
Der(3)= (Value3cont([x(1),x(2),x(3)+h])-Value3cont([x(1),x(2),x(3)-h]))/(2*h);
DiffDerivatives=Der-BelObjectiveUncondGradNAGBGP(3,x,1)';


end