clear all

NumSim=5000;
rHist0=rand(NumSim,1);
btild_1=-2;
s_=1;
InitialConditions.s0=s_;


InitialConditions.s0=1;
load ~/Golosov-Sargent/Data/Draft/cTFP.mat
[ InitialConditions.x0,InitialConditions.R0] = solveTime0Problem(Para,c,V,btild_1,s_);
[SimData(1)]=RunSimulationsFromT1AltThetaShocksUsingPolicyRules(PolicyFunctions,InitialConditions,NumSim,Para,rHist0);


load ~/Golosov-Sargent/Data/Draft/cTFPIneq.mat
[ InitialConditions.x0,InitialConditions.R0] = solveTime0Problem(Para,c,V,btild_1,s_);
[SimData(2)]=RunSimulationsFromT1AltThetaShocksUsingPolicyRules(PolicyFunctions,InitialConditions,NumSim,Para,rHist0);


load ~/Golosov-Sargent/Data/Draft/cTFPIneqBetaShocks.mat
[ InitialConditions.x0,InitialConditions.R0] = solveTime0Problem(Para,c,V,btild_1,s_);
[SimData(3)]=RunSimulationsFromT1AltThetaShocksUsingPolicyRules(PolicyFunctions,InitialConditions,NumSim,Para,rHist0);


load ~/Golosov-Sargent/Data/Draft/cTFPIneqLargerBetaShocks.mat
[ InitialConditions.x0,InitialConditions.R0] = solveTime0Problem(Para,c,V,btild_1,s_);
[SimData(4)]=RunSimulationsFromT1AltThetaShocksUsingPolicyRules(PolicyFunctions,InitialConditions,NumSim,Para,rHist0);
save ('~/Golosov-Sargent/Data/Draft/Simulations.mat','SimData','rHist0')
break;

a=2;
b=2;
figure()
subplot(2,1,1)
plot(-SimData.btildHist(a:end-b)./simData.YHist(a:end-b))
xlabel('time')
ylabel('Debt-Gdp')
subplot(2,1,2)
plot(SimData.TauHist(a:end-b))
xlabel('time')
ylabel('Labor Tax rate')



figure()
subplot(2,1,1)
plot([SimData.xHist])
xlabel('time')
ylabel('x_t')
title('Marginal utility adjusted difference in debt')
legend('IID','Persistent','Beta Shocks')
subplot(2,1,2)
plot([SimData.RHist])
xlabel('time')
ylabel('\rho_t')
title('Ratio of marginal utility')