function [ yHist sHist shockHist] = SimulateFromT1Linear(y0,G1,impact,shocks,P,T,rHist)
%SIMULATEFROMT1LINEAR Summary of this function goes here
%   Detailed explanation goes here
    yHist = zeros(length(y0),T);
    shockHist = zeros(3,T);
    sHist = zeros(1,T);
    %Assuming iid for now
    cumP = cumsum(P(:,1));
    
    sHist(1) = sum(~(rHist(1) < cumP))+1;
    yHist(:,1) = G1*y0+impact*shocks(:,sHist(1));
    shockHist(:,1) = shocks(:,sHist(1));
    
    for t= 2:T
        sHist(t) = sum(~(rHist(t) < cumP))+1;
        yHist(:,t) = G1*yHist(:,t-1)+impact*shocks(:,sHist(t));
        shockHist(:,t) = shocks(:,sHist(t));
    end
    

end

