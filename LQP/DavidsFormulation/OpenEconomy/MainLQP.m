clear all
% Params
z1bar = 3.3;
z2bar = 1;
gamma = 2;
sigma = 1;
gbar = 0.3;
bbar = -1;
beta = .9;
alpha1 = .69;
alpha2 = 1-alpha1;

% Params
Param.ss_theta_1=z1bar;
Param.ss_theta_2=z2bar;
Param.ss_g=gbar;
Param.sigma=sigma;
Param.gamma=gamma;
Param.alpha_1=alpha1;
Param.alpha_2=alpha2;
Param.beta=beta;

% Initial b2,rhof
ss_rho=3;
ss_b2=-1;

% Allocation
x0=[ss_rho^(sigma) 1 .5 .5 -.5 .5 -.05 -.05];
x=fsolve(@(x) ComputeSteadyStateAllocation_Multipliers(x,ss_b2,ss_rho,Param) ,x0);
SSVar.ss_c1=x(1);
SSVar.ss_c2=x(2);
SSVar.ss_l1=x(3);
SSVar.ss_l2=x(4);
SSVar.ss_b2=ss_b2;
SSLambda.ss_lambda_I=x(5);
SSLambda.ss_lambda_W=x(6);
SSLambda.ss_lambda_B1=x(7);
SSLambda.ss_lambda_B2=x(8);
SSLambda.ss_lambda_R=0;

% MAtrices

[Qtemp,B]=ComputeMatrixQB(SSVar,SSLambda,Param);
Q=Qtemp(1:end-1,1:end-1);
B=B(1:end-1,:);
printmat(Q,'Q','c1 c2 l1 l2 b2','c1 c2 l1 l2 b2')
printmat(B,'B','c1 c2 l1 l2 b2','g theta1 theta2')
save('Q')
save('SSVar')
save('SSLambda')


Qc1c1=Q(1,1);
Qc1l1=Q(1,3);
Qc1z1=B(1,2);
Ql1l1=Q(3,2);
Ql1z1=B(3,2);

Qc2c2=Q(2,2);
Qc2l2=Q(2,4);
Qc2z2=B(2,3);
Ql2l2=Q(4,4);
Ql2z2=B(4,3);


c1bar=SSVar.ss_c1;
c2bar=SSVar.ss_c2;
l1bar=SSVar.ss_l1;
l2bar=SSVar.ss_l2;
b2bar=SSVar.ss_b2;

taubar=1-c1bar^(sigma)*l1bar^(gamma)/z1bar
Tbar=c2bar-z2bar*(1-taubar)*l2bar
ybar=c1bar+c2bar+gbar
phi1 = l1bar^(gamma)*c1bar^sigma;
phi2 = l2bar^gamma*c2bar^sigma;
rhobar = c2bar^(-sigma)/c1bar^(-sigma);
I1 = phi1*l1bar; %FIX THIS
I2 = phi2*l2bar; %FIX This
