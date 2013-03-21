
function res=ComputeSteadyStateAllocation(x,ss_b2,ss_rho,Param)
ss_theta_1=Param.ss_theta_1;
ss_theta_2=Param.ss_theta_2;
ss_g=Param.ss_g;
sigma=Param.sigma;
gamma=Param.gamma;
beta=Param.beta;
ss_Q=beta;
ss_T=x(1);
ss_tau=x(2);
ss_c1=x(3);
ss_c2=x(4);
SSVar;
res=vertcat(res,(ss_c2^(-sigma)-ss_rho*ss_c1^(-sigma)));

end
