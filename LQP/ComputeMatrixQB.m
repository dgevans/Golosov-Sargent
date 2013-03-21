

function [Q,B]=ComputeMatrixQB(SSVar,SSLambda,Param)
ss_c1=SSVar.ss_c1;
ss_c2=SSVar.ss_c2;
ss_l1=SSVar.ss_l1;
ss_l2=SSVar.ss_l2;
ss_Q=SSVar.ss_Q;
ss_b2=SSVar.ss_b2;

ss_lambda_I=SSLambda.ss_lambda_I;
ss_lambda_R=SSLambda.ss_lambda_R;
ss_lambda_W=SSLambda.ss_lambda_W;
ss_lambda_B1=SSLambda.ss_lambda_B1;
ss_lambda_B2=SSLambda.ss_lambda_B2;

ss_theta_1=Param.ss_theta_1;
ss_theta_2=Param.ss_theta_2;
sigma=Param.sigma;
gamma=Param.gamma;
alpha_1=Param.alpha_1;
alpha_2=Param.alpha_2;
beta=Param.beta;
phi_1=ss_c1^(sigma)*ss_l1^(gamma);
phi_2=ss_c2^(sigma)*ss_l2^(gamma);

Q(:,1)=vertcat(1/2*phi_1*sigma^2*ss_lambda_W*ss_theta_2 +1/2*beta*sigma^2*ss_Q*ss_lambda_B1/ss_c1^sigma -1/2*beta*sigma^2*ss_lambda_B1/ss_c1^sigma +1/2*(ss_c1^sigma*sigma^2*ss_l1^(gamma + 1) - ss_c1)*ss_lambda_I +1/2*ss_c1*ss_lambda_R - 1/2*(ss_c2^sigma*alpha_1*sigma -ss_c2^sigma*alpha_1)*ss_c1/(ss_c1^sigma*ss_c2^sigma), 0, 1/2*(gamma +1)*ss_c1^sigma*sigma*ss_lambda_I*ss_l1^(gamma + 1) +1/2*gamma*phi_1*sigma*ss_lambda_W*ss_theta_2, 0, 0,-1/2*beta*sigma*ss_Q*ss_lambda_B1/ss_c1^sigma);

Q(:,2)=vertcat(0, -1/2*phi_2*sigma^2*ss_lambda_W*ss_theta_1 +1/2*beta*sigma^2*ss_Q*ss_lambda_B2/ss_c1^sigma -1/2*beta*sigma^2*ss_lambda_B2/ss_c2^sigma -1/2*(ss_c2^sigma*sigma^2*ss_l2^(gamma + 1) - ss_c2)*ss_lambda_I +1/2*ss_c2*ss_lambda_R - 1/2*(ss_c1^sigma*alpha_2*sigma -ss_c1^sigma*alpha_2)*ss_c2/(ss_c1^sigma*ss_c2^sigma), 0, -1/2*(gamma +1)*ss_c2^sigma*sigma*ss_lambda_I*ss_l2^(gamma + 1) -1/2*gamma*phi_2*sigma*ss_lambda_W*ss_theta_1, 0,-1/2*beta*sigma*ss_Q*ss_lambda_B2/ss_c1^sigma);

Q(:,3)=vertcat(1/2*(gamma + 1)*ss_c1^sigma*sigma*ss_lambda_I*ss_l1^(gamma + 1) +1/2*gamma*phi_1*sigma*ss_lambda_W*ss_theta_2, 0,1/2*gamma^2*phi_1*ss_lambda_W*ss_theta_2 + 1/2*(gamma^2 + 2*gamma +1)*ss_c1^sigma*ss_lambda_I*ss_l1^(gamma + 1) -1/2*ss_l1*ss_lambda_R*ss_theta_1 +1/2*(ss_c1^sigma*ss_c2^sigma*alpha_1*gamma*ss_l1^(gamma + 1) +ss_c1^sigma*ss_c2^sigma*alpha_1*ss_l1^(gamma +1))/(ss_c1^sigma*ss_c2^sigma), 0, 0, 0);


Q(:,4)=vertcat(0, -1/2*(gamma + 1)*ss_c2^sigma*sigma*ss_lambda_I*ss_l2^(gamma + 1) -1/2*gamma*phi_2*sigma*ss_lambda_W*ss_theta_1, 0,-1/2*gamma^2*phi_2*ss_lambda_W*ss_theta_1 - 1/2*(gamma^2 + 2*gamma +1)*ss_c2^sigma*ss_lambda_I*ss_l2^(gamma + 1) -1/2*ss_l2*ss_lambda_R*ss_theta_2 +1/2*(ss_c1^sigma*ss_c2^sigma*alpha_2*gamma*ss_l2^(gamma + 1) +ss_c1^sigma*ss_c2^sigma*alpha_2*ss_l2^(gamma +1))/(ss_c1^sigma*ss_c2^sigma), 0, 0);

Q(:,5)=vertcat(0, 0, 0, 0, 0, 1/2*beta*ss_b2*ss_lambda_I);

Q(:,6)=vertcat(-1/2*beta*sigma*ss_Q*ss_lambda_B1/ss_c1^sigma,-1/2*beta*sigma*ss_Q*ss_lambda_B2/ss_c1^sigma, 0, 0,1/2*beta*ss_b2*ss_lambda_I, 1/2*beta*ss_b2*ss_lambda_I +1/2*(ss_Q*ss_lambda_B1/ss_c1^sigma +ss_Q*ss_lambda_B2/ss_c1^sigma)*beta);



B(:,1)=vertcat(0, 0, 0, 0, 0, 0);
B(:,2)=vertcat(0, -1/2*phi_2*sigma*ss_lambda_W*ss_theta_1,-1/2*ss_l1*ss_lambda_R*ss_theta_1,-1/2*gamma*phi_2*ss_lambda_W*ss_theta_1, 0, 0);
B(:,3)=vertcat(1/2*phi_1*sigma*ss_lambda_W*ss_theta_2, 0,1/2*gamma*phi_1*ss_lambda_W*ss_theta_2,-1/2*ss_l2*ss_lambda_R*ss_theta_2, 0, 0)  ;
end

