
    function Tau=GetTaxRates(x,R,s_,c,V,domain,PolicyRulesStore,Para,S)
        [PolicyRulesInit]=GetInitialApproxPolicy([x R s_] ,domain,PolicyRulesStore);
        [PolicyRules, ~,~,~]=CheckGradNAG(x,R,s_,c,V,PolicyRulesInit,Para);
     c1=PolicyRules(1:S);
     uc1=Para.psi./(c1.^(Para.sigma));
    l1=PolicyRules(2*S+1:3*S);
    ul1=(1-Para.psi)./(1-l1);
    Tau=1-(ul1./(Para.theta_1.*uc1));
    end
    