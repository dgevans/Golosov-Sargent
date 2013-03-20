function [ x,R,PolicyRule ] = findSteadyStateN( x0,R0,Para)
%FINDSTEADYSTATE Summary of this function goes here
%   Detailed explanation goes here
      S = length(Para.P);
      N = Para.N;
      options = optimset('Display','off','TolFun',1e-14);
      nMult = (4*S+4)*N-(S+4) - ( (2*S+2)*N-2);
      X0 = [0.5*ones(2*S*N,1);x0;R0;zeros(nMult,1)];
      [PolicyRule,~,exitFlag] = fsolve(@(Xtemp)SSResidualsN(Xtemp,Para),X0,options);
      if(exitFlag <1)
          exception = MException('Failed to find Steady State');
          throw(exception);
      end
          
     
      x = PolicyRule(2*S*N+1:(2*S+1)*N-1)';
      R = PolicyRule((2*S+1)*N:(2*S+2)*N-2)';
end


function [res] = findMultipliers(X,Mult,Para)
    S = length(Para.P);
    Xtemp = [X;Mult];
    nMult = (4*S+4)*N-(2*S+4) - ( (2*S+2)*N-2);
    
    resTemp = SSResidualsN(Xtemp,Para);
    res(1:nMult,1) = resTemp((2*S+2)*N-1:9*S);
end


function [res] = SSResidualsTruncated(Xtild,R0,x0,Para)
    N = Para.N;
    S = length(Para.P);
    
    nMult = (4*S+4)*N-(S+4) - ( (2*S+2)*N-2);
    
    X = [Xtild;R0;x0;zeros(nMult,1)];
    
    res = SSResidualsN(X,Para);
    res = res(1:2*S*N);
end