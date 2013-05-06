
function res=ComputeSteadyStateAllocation_Multipliers(x,ss_b,Param)
ss_g=Param.ss_g;
sigma=Param.sigma;
gamma=Param.gamma;
beta=Param.beta;
theta=Param.theta;
ss_c=x(1);
ss_l=x(2);
ss_lambda_I=x(3);
ss_lambda_B=x(4);
ss_lambda_R=x(5);
ss_Q=beta;
SSVar;
SSLambda;
res=vertcat(res_ss,res_lambda);
end
