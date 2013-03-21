clear all
clc
% Params
Param.ss_theta_1=3.3;
Param.ss_theta_2=1;
Param.ss_g=.15;
Param.sigma=1;
Param.gamma=1;
Param.alpha_1=.69;
Param.alpha_2=.31;
Param.beta=.9;

% Initial b2,rho
ss_rho=2.8;
ss_b2=-1;

% Allocation
x0=[.2 .3 .5 .5];
x=fsolve(@(x) ComputeSteadyStateAllocation(x,ss_b2,ss_rho,Param) ,x0);
SSVar.ss_tau=x(1);
SSVar.ss_T=x(2);
SSVar.ss_c1=x(3);
SSVar.ss_c2=x(4);
SSVar.ss_Q=Param.beta;
SSVar.ss_b2=ss_b2;
SSVar
% Multipliers
x0=[.5 .5 .5 .5 .5];
x=fsolve(@(x) ComputeSSMultipliers(x,SSVar,Param) ,x0);

SSLambda.ss_lambda_I=x(1);
SSLambda.ss_lambda_R=x(2);
SSLambda.ss_lambda_B2=x(3);
SSLambda.ss_lambda_E1=x(4);
SSLambda.ss_lambda_E2=x(5);

% MAtrices

[Qtemp,B]=ComputeMatrixQB(SSVar,SSLambda,Param);
Q=Qtemp(1:end-1,1:end-1)
B=B(1:end-1,:)
printmat(Q,'Q','tau T c1 c2 b2','tau T c1 c2 b2')
printmat(B,'B','tau T c1 c2 b2','g theta1 theta2')
Qtautau=Q(1,1);
Qtauc1=Q(1,3);
Qtauc2=Q(1,4);
QTT=Q(2,2);
Qc1c1=Q(3,3);
Qc2c2=Q(4,4);
Btheta_1tau=B(1,2);
Btheta_2tau=B(1,3);
Btheta_1c1=B(3,2);
Btheta_2c2=B(4,3);

x0=[]
x0(1)=Q(3,3);
x0(2)=Q(4,4);

x=fsolve(@(x) GetCoeff(x,Q,B),x0);

omega_1=x(1)
omega_2=x(2)

eta_1=sqrt(omega_1/Qc1c1)*Qtauc1;
eta_2=sqrt(omega_2/Qc2c2)*Qtauc2;

kappa_1=sqrt(omega_1/Qc1c1)*Btheta_1c1;
kappa_2=sqrt(omega_2/Qc2c2)*Btheta_2c2;

delta_1=sqrt(Qc1c1/omega_1)
delta_2=sqrt(Qc2c2/omega_2)
omega_3=Qtautau-omega_1*eta_1^2-omega_2*eta_2^2
omega_4=QTT
