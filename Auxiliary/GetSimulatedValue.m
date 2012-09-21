function [ V ] = GetSimulatedValue(CoeffFileName,u2btild0,R0,NumSim,T,Para)
%GETSIMULATEDVALUE Summary of this function goes here
%   Detailed explanation goes here

load(CoeffFileName)

VSim = zeros(NumSim,1);

alpha = [Para.alpha_1,Para.alpha_2];
psi = Para.psi;
beta = Para.beta;

for i = 1:NumSim
    disc = 1;
    u2btild = u2btild0;
    R = R0;
    s_ = 1;
    for t = 1:T
        [PolicyRulesInit]=GetInitialApproxPolicy([u2btild R s_] ,x_state,PolicyRulesStore);
        [PolicyRules, ~,exitflag,~]=CheckGradNAG(u2btild,R,s_,c,V,PolicyRulesInit,Para,0);
        
        %Get Policies
        c1=PolicyRules(1:2);
        c2=PolicyRules(3:4);
        l1=PolicyRules(5:6);
        l2=PolicyRules(7:8);
        Rprime=PolicyRules(end-3:end-2);
        % x' - u_c_2* btildprime
        u2btildprime=PolicyRules(end-1:end);
        
        %get next state
        s = randi(2);
        
        %update value
        VSim = VSim + disc *(alpha(1)*uBGP(c1(s),l1(s),psi)+alpha(2)*uBGP(c2(s),l2(s),psi));
        
        %Store next period states
        u2btild = u2btildprime(s);
        R = Rprime(s);
        s_ = s;
        disc = disc *beta;
        
    end
    i
end

V = mean(VSim);
end

