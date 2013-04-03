clear all
clc
% Params
z1bar = 3.3;
z2bar = 1;
gamma = 1;
sigma = 1;
gbar = 0.3484;
bbar = -1.0471;
beta = .9;
psi=0.6958;
alpha1 = .69;
alpha2 = 1-alpha1;
% Initial b2,rhof
ss_rho=3.1433;
ss_b2=bbar;

% Params
Param.psi=psi;
Param.ss_theta_1=z1bar;
Param.ss_theta_2=z2bar;
Param.ss_g=gbar;
Param.sigma=sigma;
Param.gamma=gamma;
Param.alpha_1=alpha1;
Param.alpha_2=alpha2;
Param.beta=beta;
Param.bgp=1
% Allocation
x0=[ss_rho 1 .7 .5];
x=fsolve(@(x) ComputeSteadyStateAllocation(x,ss_b2,ss_rho,Param) ,x0);
SSVar.ss_c1=x(1);
SSVar.ss_c2=x(2);
SSVar.ss_l1=x(3);
SSVar.ss_l2=x(4);
SSVar.ss_Q=Param.beta;
SSVar.ss_b2=ss_b2;

% Multipliers
x0=[-.5 -.5 .5 .5 .5];
options=optimset('TolX',1e-8);
x=fsolve(@(x) ComputeSSMultipliers(x,SSVar,Param) ,x0,options);

SSLambda.ss_lambda_I=x(1);
SSLambda.ss_lambda_R=x(2);
SSLambda.ss_lambda_W=x(3);
SSLambda.ss_lambda_B1=x(4);
SSLambda.ss_lambda_B2=x(5);

% MAtrices
format short
[Omega,B]=ComputeMatrixQB(SSVar,SSLambda,Param);
printmat(Omega,'Omega','c1 c2 l1 l2 b2 Q','c1 c2 l1 l2 b2 Q')
printmat(B,'B','c1 c2 l1 l2 b2 Q','g theta1 theta2')


Omegac1c1=Omega(1,1);
Omegac1l1=Omega(1,3);
Omegac1z1=B(1,2);
Omegal1l1=Omega(3,3);
Omegal1z1=B(3,2);
Omegac1Q=Omega(1,6);
Omegac2c2=Omega(2,2);
Omegac2l2=Omega(2,4);
Omegac2z2=B(2,3);
Omegal2l2=Omega(4,4);
Omegal2z2=B(4,3);
Omegac2Q=Omega(2,6);
Omegab2Q=Omega(5,6);

c1bar=SSVar.ss_c1;
c2bar=SSVar.ss_c2;
l1bar=SSVar.ss_l1;
l2bar=SSVar.ss_l2;
b2bar=SSVar.ss_b2;

uc1=psi*c1bar^(-sigma);
ul1=(1-psi)/(1-l1bar);
uc2=psi*c2bar^(-sigma);
ul2=(1-psi)/(1-l2bar);




taubar=1-ul1/(z1bar*uc1);
Tbar=c2bar-z2bar*(1-taubar)*l2bar;
ybar=c1bar+c2bar+gbar;
phi1 = ul1/uc1;
phi2 = ul2/uc2;
rhobar = c2bar^(-sigma)/c1bar^(-sigma);
I1 = phi1*l1bar; %FIX THIS
I2 = phi2*l2bar; %FIX This
gamma1_I = 1/(1-l1bar);
gamma2_I = 1/(1-l2bar);
gamma1 = l1bar*gamma1_I;
gamma2 = l2bar*gamma2_I;