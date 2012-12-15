% This script checks the necessary condition for existence of steady state
clear all
clc
options = optimset('Display','off','TolX',1e-10);
                
load('Data/temp/csigmalow.mat')
 xState=fsolve(@(x) GetCrossingPoints(x,1,c,V,PolicyRulesStore,x_state,Para),[0 4],options)
u2btild_=xState(1);
R_=xState(2);

s_=1;
[xInit]=GetInitialApproxPolicy([u2btild_ R_ s_],x_state,PolicyRulesStore);
[PolicyRules, V_new,exitflag,~]=CheckGradNAG(u2btild_,R_,s_,c,V,xInit',Para,0);
BelObjectiveUncondGradNAGBGP(3,PolicyRules(1:3),1,1)
xprime=[u2btild_ u2btild_];
Rprime=[R_ R_]
    %compute objective
    V_x(:,1)=funeval(c(s_,:)',V(1),[xprime(1,1) Rprime(1,1)],[1,0]);
    V_x(:,2)=funeval(c(s_,:)',V(2),[xprime(1,2) Rprime(1,2)],[1,0]);
    V_R(:,1)=funeval(c(s_,:)',V(1),[xprime(1,1) Rprime(1,1)],[0,1])
    V_R(:,2)=funeval(c(s_,:)',V(2),[xprime(1,2) Rprime(1,2)],[0,1])




c2xR=fsolve(@(c2xR) ResNecessaryConditionForSS(c2xR,s_,Para) ,[PolicyRules(3:4) xState],options)

                cRat = R_^(-1/Para.sigma);
                c1_1 = (0.8*(Para.n1*Para.theta_1+Para.n2*Para.theta_2)-Para.g(1))/(Para.n1+cRat*Para.n2);
                c1_2 = (0.8*(Para.n1*Para.theta_1+Para.n2*Para.theta_2)-Para.g(2))/(Para.n1+cRat*Para.n2);
                c2_1 = cRat*c1_1;
                [xSS,~,exitFlag] = fsolve(@(x) SteadyStateResiduals(x,u2btild_,R_,Para,s_),[c1_1 c1_2 c2_1],options);
                [res, c1_, c2_, l1_, l2_] = SteadyStateResiduals(xSS,u2btild_,R_,Para,s_);

ResNecessaryConditionForSS([c2_ xState],s_,Para)
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
