% Inputs - xInit, state variables - u2btild,,R,s_  coeff, value
% function, para
function [PolicyRules, V_new,exitflag,fvec]=CheckGradNAG(u2bdiff,RR,s,c,VV,xInit,Para,flagOpt)
global V Vcoef R u2btild Par s_ flagCons
Para.epsilon=1e-5;
Para.k_a=2;
%Get the initial guess for the uconstraint problem. With the simplification
%we need only c1_1,c1_2and c2_1

xInit=xInit(1:3);
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

%% Now solve the unconstraint problem FOC using NAG
% use the last solution
warning('off', 'NAG:warning')
[x, fvec,exitflag]=c05nb('BelObjectiveUncondGradNAGBGP',xInit);
if exitflag==4
    exitflag=-2;
else
    exitflag=1;
end


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

% State next period
X(1,:) = [psi*c2_1^(-1)*btildprime(1),c2_1^(-1)/c1_1^(-1)];%state next period
X(2,:) = [psi*c2_2^(-1)*btildprime(2),c2_2^(-1)/c1_2^(-1)];%state next period


%compute objective
if ~isreal(X)
    X=real(X);
end
Vobj = P(s_,1)*(alpha(1)*uBGP(c1_1,l1(1),psi)+alpha(2)*uBGP(c2_1,l2(1),psi)...
    +beta*funeval(Vcoef{1},V(1),X(1,:)));

Vobj = Vobj + P(s_,2)*(alpha(1)*uBGP(c1_2,l1(2),psi)+alpha(2)*uBGP(c2_2,l2(2),psi)...
    +beta*funeval(Vcoef{2},V(2),X(2,:)));

V_new=Vobj;
PolicyRules=[c1_1 c1_2 c2_1 c2_2 l1(1) l1(2) l2(1) l2(2) btildprime c2_1^(-1)/c1_1^(-1) c2_2^(-1)/c1_2^(-1) u2btildprime(1) u2btildprime(2)];
end
