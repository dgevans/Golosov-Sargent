% This script checks the necessary condition for existence of steady state
clear all

load('Data/temp/csigmamed.mat')
 xState=fsolve(@(x) GetCrossingPoints(x,1,c,V,PolicyRulesStore,x_state,Para),[.8 4])
u2btild=xState(1);
R=xState(2);
s_=1;
[xInit]=GetInitialApproxPolicy([u2btild R s_],x_state,PolicyRulesStore);
[PolicyRules, V_new,exitflag,~]=CheckGradNAG(u2btild,R,s_,c,V,xInit',Para,0);

fsolve(@(c2xR) ResNecessaryConditionForSS(c2xR,s_,Para) ,[PolicyRules(3:4) xState])
% 
%   psi = Para.psi;
%     sigma = Para.sigma;
%     beta =  Para.beta;
%     P = Para.P;
%     theta_1 = Para.theta_1;
%     theta_2 = Para.theta_2;
%     g = Para.g;
%     alpha_2 = Para.alpha_2;
%     n1 = Para.n1;
%     n2 = Para.n2;
% c1=PolicyRules(1:2); 
% c2=PolicyRules(3:4);
% l1=PolicyRules(5:6);
% l2=PolicyRules(5:6);
% uc2=psi*c2.^(-sigma);
% ucc2=-sigma*psi*c2.^(-sigma-1);
%  Euc2 = psi*c2.^(-sigma)*(P(s_,:)');   
% Term1=alpha_2*c2.^(-sigma);
% Term2=(1-R^(1/sigma))*psi*(1-sigma)*c2.^(-sigma);
% Term3=((1-psi)*R./(1-l1).^2)*(1/(theta_2*(n2+R*n1)))*(R*theta_2/theta_1);
% Term4=((1-psi)./((1-l2).^2)).*(1/(theta_2*(n2+R*n1)));
% Term5=u2btild/beta;
% Term6=uc2.*P(s_,:).*ucc2*(-1/(Euc2)^2)+ (1/Euc2)*ucc2;
% Term1./(Term2+(Term3-Term4)-Term5.*Term6)
