function res=ComputeSSMultipliers(x,SSVar,Param)

ss_lambda_I=x(1);
ss_lambda_R=x(2);
ss_lambda_B2=x(3);
ss_lambda_E1=x(4);
ss_lambda_E2=x(5);

ss_tau=SSVar.ss_tau;
ss_T=SSVar.ss_T;
ss_c1=SSVar.ss_c1;
ss_c2=SSVar.ss_c2;
ss_Q=SSVar.ss_Q;
ss_b2=SSVar.ss_b2;



ss_theta_1=Param.ss_theta_1;
ss_theta_2=Param.ss_theta_2;
ss_g=Param.ss_g;
sigma=Param.sigma;
gamma=Param.gamma;
alpha_1=Param.alpha_1;
alpha_2=Param.alpha_2;
beta=Param.beta;

SSLambda;
end



