function [ x,R,PolicyRule ] = findSteadyState( x0,R0,Para)
%FINDSTEADYSTATE Summary of this function goes here
%   Detailed explanation goes here
      S = length(Para.P);
      if min(Para.theta_2) > 0
            cRat = R0^(-1/Para.sigma);
            c1 = (0.8*(Para.n1*Para.theta_1+Para.n2*Para.theta_2)-Para.g)/(Para.n1+cRat*Para.n2);
            c2_ = cRat*c1; c2_(S) = [];

            options = optimset('Display','off','TolFun',1e-10);
            [xSS,~,~] = fsolve(@(x) SteadyStateResiduals(x,x0,R0,Para,1),[c1 c2_],options);
            [~, c1_, c2_, l1_, l2_] = SteadyStateResiduals(xSS,x0,R0,Para,1);

            X0 = [c1_,c2_,l1_,l2_,x0,R0];
            nMult = 3*S+2;

            f = @(Mult) findMultipliers(X0,Mult,Para);
            I = eye(nMult);

            b = -f(zeros(1,nMult))';
            for i = 1:nMult
                A(:,i) = f(I(i,:))'+b;
            end
            Mult = A\b;
            X0(4*S+3:7*S+4) = Mult;
      else
          nMult = 2*S+2;
          X0 = [0.5*ones(3*S,1);x0;R0;zeros(nMult,1)]';
      end
      options = optimset('Display','off','TolFun',1e-14);
      [PolicyRule,~,exitFlag] = fsolve(@(Xtemp)SSResiduals(Xtemp,Para),X0,options);
      if(exitFlag <1)
          exception = MException('Failed to find Steady State');
          throw(exception);
      end
          
      if(min(Para.theta_2) > 0)
        R = PolicyRule(4*S+2);
        x = PolicyRule(4*S+1);
      else
        x = PolicyRule(3*S+1);
        R = PolicyRule(3*S+2);
      end
      
end

function [fvec, user, iflag] = SSResidualsNAG(n,x,user,iflag)
    fvec = SSResiduals(x,user.Para);
end

function [res] = findMultipliers(X,Mult,Para)
    S = length(Para.P);
    nMult = 3*S+2;
    X(4*S+3:7*S+4) = Mult;
    
    resTemp = SSResiduals(X,Para);
    res(1:nMult) = resTemp(6*S-1:9*S);
end
