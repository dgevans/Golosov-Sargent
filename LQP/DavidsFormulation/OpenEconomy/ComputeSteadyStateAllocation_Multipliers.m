
function res=ComputeSteadyStateAllocation_Multipliers(x,ss_b2,ss_rho,Param)
ss_theta_1=Param.ss_theta_1;
ss_theta_2=Param.ss_theta_2;
ss_g=Param.ss_g;
sigma=Param.sigma;
gamma=Param.gamma;
beta=Param.beta;
alpha_1=Param.alpha_1;
alpha_2=Param.alpha_2;
ss_c1=x(1);
ss_c2=x(2);
ss_l1=x(3);
ss_l2=x(4);
ss_lambda_I=x(5);
ss_lambda_W=x(6);
ss_lambda_B1=x(7);
ss_lambda_B2=x(8);
ss_lambda_R=0;
phi_1=ss_c1^(sigma)*ss_l1^(gamma);
phi_2=ss_c2^(sigma)*ss_l2^(gamma);
SSVar;
SSLambda;
res_rho=ss_c2^(-sigma)-ss_rho*ss_c1^(-sigma);
res=vertcat(res_ss,res_lambda,res_rho);
end
