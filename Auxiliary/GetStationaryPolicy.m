function [ PolicyRules] = GetStationaryPolicy(xState,Para)
s=1;
u2bdiff=xState(1);
RR=xState(2);
options=optimset('Display','off');
[x,fvec,exitval]=fsolve(@(x) SteadyStateCharacterization(x,u2bdiff,RR,Para,1) ,[1 1 RR^-1],options);
Para.theta=[Para.theta_1 Para.theta_2];
Para.alpha=[Para.alpha_1 Para.alpha_2];
Par=Para;
u2btild=u2bdiff;
R=RR;
s_=s;
n1=Para.n1;
n2=Para.n2;
ctol=Para.ctol;

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
Vobj = P(s_,1)*(alpha(1)*uBGP(c1_1,l1(1),psi)+alpha(2)*uBGP(c2_1,l2(1),psi));
Vobj = Vobj + P(s_,2)*(alpha(1)*uBGP(c1_2,l1(2),psi)+alpha(2)*uBGP(c2_2,l2(2),psi));
PolicyRules=[c1_1 c1_2 c2_1];


end

