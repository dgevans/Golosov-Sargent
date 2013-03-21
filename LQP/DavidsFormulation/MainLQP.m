clear all
% Params
Param.ss_theta_1=3.3;
Param.ss_theta_2=1;
Param.ss_g=.3;
Param.sigma=1;
Param.gamma=1;
Param.alpha_1=.69;
Param.alpha_2=.31;
Param.beta=.9;

% Initial b2,rhof
ss_rho=3;
ss_b2=1;

% Allocation
x0=[.5 .5 .5 .5];
x=fsolve(@(x) ComputeSteadyStateAllocation(x,ss_b2,ss_rho,Param) ,x0);
SSVar.ss_c1=x(1);
SSVar.ss_c2=x(2);
SSVar.ss_l1=x(3);
SSVar.ss_l2=x(4);
SSVar.ss_Q=Param.beta;
SSVar.ss_b2=ss_b2;

% Multipliers
x0=[.5 .5 .5 .5 .5];
x=fsolve(@(x) ComputeSSMultipliers(x,SSVar,Param) ,x0);

SSLambda.ss_lambda_I=x(1);
SSLambda.ss_lambda_R=x(2);
SSLambda.ss_lambda_W=x(3);
SSLambda.ss_lambda_B1=x(4);
SSLambda.ss_lambda_B2=x(5);

% MAtrices

[Qtemp,B]=ComputeMatrixQB(SSVar,SSLambda,Param);
Q=Qtemp(1:end-1,1:end-1);
B=B(1:end-1,:);
printmat(Q,'Q','c1 c2 l1 l2 b2','c1 c2 l1 l2 b2')
printmat(B,'B','c1 c2 l1 l2 b2','g theta1 theta2')
