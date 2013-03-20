% Inputs - xInit, state variables - u2btild,,R,s_  coeff, value
% function, para
function [res]=StationaryFOC(StationaryXR,s,c,VV,Para)
global V Vcoef R u2btild Par s_ 
u2bdiff=StationaryXR(1);
RR=StationaryXR(2);
options=optimset('Display','off');
[x,fvec,exitval]=fsolve(@(x) SteadyStateCharacterization(x,u2bdiff,RR,Para,s) ,[1 1 RR^-1],options);
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

res(1:3)=BelObjectiveUncondGradNAGBGP(3,x,1);

%% GET THE Policy Rules
psi= Par.psi;
beta =  Par.beta;
P = Par.P;
theta_1 = Par.theta(1);
theta_2 = Par.theta(2);
g = Par.g;
alpha = Par.alpha;

sigma = 1;
c1_1=x(1);
c1_2=x(2);
c2_1=x(3);

%compute components from unconstrained guess
[c2_2 grad_c2_2] = computeC2_2(c1_1,c1_2,c2_1,R,s_,P,sigma);
[l1 l1grad l2 l2grad] = computeL(c1_1,c1_2,c2_1,c2_2,grad_c2_2,...
    theta_1,theta_2,g,n1,n2);
[btildprime grad_btildprime] = computeBtildeprime(c1_1,c1_2,c2_1,c2_2,grad_c2_2,l1,l2,l1grad,l2grad,...
    u2btild,s_,psi,beta,P);

% x' - definition
u2btildprime=psi*[c2_1^(-1) c2_2^(-1)].*btildprime;
Rprime(1)=c1_1./c2_1;
Rprime(2)=c1_2./c2_2;
end
