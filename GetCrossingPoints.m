% Inputs - xInit, state variables - u2btild,,R,s_  coeff, value
% function, para
function [res]=GetCrossingPoints(x,s_,c,V,PolicyRulesStore,x_state,Para)
n1=Para.n1;
n2=Para.n2;
alpha_1=Para.alpha_1;
alpha_2=Para.alpha_2;
theta_1=Para.theta_1;
theta_2=Para.theta_2;
psi=Para.psi;
beta=Para.beta;
u2btild=x(1);
R=x(2);
    [PolicyRulesInit]=GetInitialApproxPolicy([u2btild R s_],x_state,PolicyRulesStore);
    [PolicyRules, V_new,exitflag,~]=CheckGradNAG(u2btild,R,s_,c,V,PolicyRulesInit,Para,0) ;
        u2BtildePrime=PolicyRules(end-1:end);
        Rprime=PolicyRules(end-3:end-2);
    res(1)=u2BtildePrime(1)-u2BtildePrime(2);
res(2)=Rprime(1)-Rprime(2);
res(3)=Rprime(1)-R;
res(2)=u2BtildePrime(1)-u2btild;

end
