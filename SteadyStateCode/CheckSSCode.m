%Get Steady state from solution
function [PolicyRules,res] = CheckSSCode()

load Data/Calibration/csigmaMed.mat
X = fsolve(@(x)GetCrossingPoints(x,1,c,V,PolicyRulesStore,x_state,Para),[1.5,4]);
u2btild = X(1);
R = X(2);
s_ = 1;
[PolicyRulesInit]=GetInitialApproxPolicy([u2btild R s_],x_state,PolicyRulesStore);
[PolicyRules, V_new,exitflag,~]=CheckGradNAG(u2btild,R,s_,c,V,PolicyRulesInit,Para,0) ;
PolicyRules(14) = [];
PolicyRules(12) = [];
PolicyRules(9:10) = [];


PolicyRules(11:18) = zeros(1,8);

PolicyRules = fsolve(@(x)SSResiduals(x,Para),PolicyRules);



end


function [res] = findMultipliers(PolicyRules,X,Para)
    x = PolicyRules;
    x(11:18) = X;
    
    resTemp = SSResiduals(x,Para);
    res(1:8) = resTemp(11:18);
end
