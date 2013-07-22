clear all

NumSim=25000;
rHist0=rand(NumSim,1);
btild_1=-2;
s_=1;
InitialConditions.s0=s_;


% Benchmark : TFP+INEQ+BETA(Persistence)
InitialConditions.s0=1;
load ~/Golosov-Sargent/Data/Draft/cTFPIneqBeta.mat
[ InitialConditions.x0,InitialConditions.R0] = solveTime0Problem(Para,c,V,btild_1,s_);
PolicyFunctions(1)=PolicyFunctions;
%[PolicyFunctions(1),errorPolicyFunctions]=GetPolicyRuleApproximations(Para,c,V,35,35,PolicyRulesStore,'spli',25,25,3,domain);
[SimData(1)]=RunSimulationsFromT1AltThetaShocksUsingPolicyRules(PolicyFunctions(1),InitialConditions,NumSim,Para,rHist0);



% TFP+INEQ

load ~/Golosov-Sargent/Data/Draft/cTFPIneq.mat
[ InitialConditions.x0,InitialConditions.R0] = solveTime0Problem(Para,c,V,btild_1,s_);
PolicyFunctions(2)=PolicyFunctions;
%[PolicyFunctions(2),errorPolicyFunctions]=GetPolicyRuleApproximations(Para,c,V,35,35,PolicyRulesStore,'spli',25,25,3,domain)
[SimData(2)]=RunSimulationsFromT1AltThetaShocksUsingPolicyRules(PolicyFunctions(2),InitialConditions,NumSim,Para,rHist0);


% TFP
load ~/Golosov-Sargent/Data/Draft/cTFP.mat
[ InitialConditions.x0,InitialConditions.R0] = solveTime0Problem(Para,c,V,btild_1,s_);
PolicyFunctions(3)=PolicyFunctions;
%[PolicyFunctions(3),errorPolicyFunctions]=GetPolicyRuleApproximations(Para,c,V,35,35,PolicyRulesStore,'spli',25,25,3,domain)
[SimData(3)]=RunSimulationsFromT1AltThetaShocksUsingPolicyRules(PolicyFunctions(3),InitialConditions,NumSim,Para,rHist0);


% Large gShcoks
load ~/Golosov-Sargent/Data/Draft/cGShocks.mat
[ InitialConditions.x0,InitialConditions.R0] = solveTime0Problem(Para,c,V,btild_1,s_);
PolicyFunctions(4)=PolicyFunctions;
%[PolicyFunctions(4),errorPolicyFunctions]=GetPolicyRuleApproximations(Para,c,V,35,35,PolicyRulesStore,'spli',25,25,3,domain)
[SimData(4)]=RunSimulationsFromT1AltThetaShocksUsingPolicyRules(PolicyFunctions(4),InitialConditions,NumSim,Para,rHist0);
save('~/Golosov-Sargent/Data/Draft/Simulations.mat','SimData','rHist0','InitialConditions','NumSim')
