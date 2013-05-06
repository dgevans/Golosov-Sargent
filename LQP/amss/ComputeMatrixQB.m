

function [Q,B]=ComputeMatrixQB(SSVar,SSLambda,Param)
ss_c=SSVar.ss_c;
ss_l=SSVar.ss_l;
beta=Param.beta;
ss_b=SSVar.ss_b;
ss_lambda_I=SSLambda.ss_lambda_I;
ss_lambda_B=SSLambda.ss_lambda_B;
ss_lambda_R=SSLambda.ss_lambda_R;
ss_Q=beta;
sigma=Param.sigma;
gamma=Param.gamma;

theta=Param.theta;
QMat;
BMat;

end

