
    function TauCor=getTaxMoments(x,R,s_,c,V,domain,PolicyRulesStore,Para,S)
        [PolicyRulesInit]=GetInitialApproxPolicy([x R s_] ,domain,PolicyRulesStore);
        [PolicyRules]=CheckGradNAG(x,R,s_,c,V,PolicyRulesInit,Para);
            
        Rprime=PolicyRules(end-2*S+1:end-S);
        % x' - u_c_2* debtprime
        xprime=PolicyRules(end-S+1:end);
        tau = GetTaxRates(x,R,s_,c,V,domain,PolicyRulesStore,Para,S)';
        tauPrime = zeros(S,S);
        w = zeros(S,S);
        for s = 1:S
            w(:,s) = Para.P(s_,s)*Para.P(s,:)';
            tauPrime(:,s) =GetTaxRates(xprime(s),Rprime(s),s,c,V,domain,PolicyRulesStore,Para,S);
        end
        tau = kron(tau,ones(S,1));
        tauPrime = reshape(tauPrime,S*S,1);
        w = reshape(w,S*S,1);
        Etau = dot(w,tau);
        EtauPrime = dot(w,tauPrime);
        varTau = dot((tau-Etau).*(tau-Etau),w);
        varTauPrime = dot((tauPrime-EtauPrime).*(tauPrime-EtauPrime),w);
        covTauTauPrime = dot((tau-Etau).*(tauPrime-EtauPrime),w);
        TauCor = covTauTauPrime/(sqrt(varTau)*sqrt(varTauPrime));
        
    end