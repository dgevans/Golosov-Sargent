function res=ComputeSSMultipliers(x,SSVar,Param)

ss_lambda_I=x(1);
ss_lambda_R=x(2);
ss_lambda_W=x(3);
ss_lambda_B1=x(4);
ss_lambda_B2=x(5);

ss_c1=SSVar.ss_c1;
ss_c2=SSVar.ss_c2;
ss_l1=SSVar.ss_l1;
ss_l2=SSVar.ss_l2;
ss_Q=SSVar.ss_Q;
ss_b2=SSVar.ss_b2;



ss_theta_1=Param.ss_theta_1;
ss_theta_2=Param.ss_theta_2;
sigma=Param.sigma;
gamma=Param.gamma;
alpha_1=Param.alpha_1;
alpha_2=Param.alpha_2;
beta=Param.beta;
phi_1=ss_c1^(sigma)*ss_l1^(gamma);
phi_2=ss_c2^(sigma)*ss_l2^(gamma);
SSLambda;

end



