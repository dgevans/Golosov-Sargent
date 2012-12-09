function [ R,x,PolicyRule ] = findSteadyState( x0,R0,Para)
%FINDSTEADYSTATE Summary of this function goes here
%   Detailed explanation goes here
      cRat = R0^(-1/Para.sigma);
      c1_1 = (0.8*(Para.n1*Para.theta_1+Para.n2*Para.theta_2)-Para.g(1))/(Para.n1+cRat*Para.n2);
      c1_2 = (0.8*(Para.n1*Para.theta_1+Para.n2*Para.theta_2)-Para.g(2))/(Para.n1+cRat*Para.n2);
      c2_1 = cRat*c1_1;
      
      options = optimset('Display','off');
      [xSS,~,~] = fsolve(@(x) SteadyStateResiduals(x,x0,R0,Para,1),[c1_1 c1_2 c2_1],options);
      [~, c1_, c2_, l1_, l2_] = SteadyStateResiduals(xSS,x0,R0,Para,1);
      
      X = [c1_,c2_,l1_,l2_,R0,x0];
      
      f = @(Mult) findMultipliers(X,Mult,Para);
      I = eye(8);
      
      b = -f(zeros(1,8))';
      for i = 1:8
          A(:,i) = f(I(i,:))'+b;
      end
      Mult = A\b;
      
      X(11:18) = Mult;
      
      
      PolicyRule = fsolve(@(Xtemp)SSResiduals(Xtemp,Para),X);
      x = PolicyRule(10);
      R = PolicyRule(9);
      
end


function [res] = findMultipliers(X,Mult,Para)
    X(11:18) = Mult;
    
    resTemp = SSResiduals(X,Para);
    res(1:8) = resTemp(11:18);
end
