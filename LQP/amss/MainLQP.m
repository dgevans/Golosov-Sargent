clc
clear all
% Params
gamma = 2;
sigma = 1;
gbar = 0.3;
bbar = -1;
beta = .9;
theta=1;
% Params
Param.theta=theta;
Param.ss_g=gbar;
Param.sigma=sigma;
Param.gamma=gamma;
Param.beta=beta;

% Initial b2,rhof
ss_b=-1;
% Allocation
ss_l0=1;
x0=[ss_l0-gbar ss_l0 -.5 -.5 -.5];
x=fsolve(@(x) ComputeSteadyStateAllocation_Multipliers(x,ss_b,Param) ,x0);
SSVar.ss_c=x(1);
SSVar.ss_l=x(2);
SSVar.ss_b=ss_b;
SSVar.ss_Q=beta;
SSLambda.ss_lambda_I=x(3);
SSLambda.ss_lambda_B=x(4);
SSLambda.ss_lambda_R=x(5);

% MAtrices

[Qtemp,B]=ComputeMatrixQB(SSVar,SSLambda,Param);
Q=Qtemp;
B=B;
printmat(Q,'Q','c l b Q','c l b Q')
printmat(B,'B','c l b Q','g')
save('Q')
save('SSVar')
save('SSLambda')


Qcc=Q(1,1);
Qcl=Q(1,3);
Qll=Q(2,2);
Qcg=B(1,1);
Qlg=B(2,1);

