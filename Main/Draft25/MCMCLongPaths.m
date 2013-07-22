clear all
matlabpool open local
load ~/Golosov-Sargent/Data/Draft/cTFPIneqLargerBetaShocks.mat
NumInit=5
NumSim=5000;
btild_1=-2;
s_=1;
InitialConditions.s0=s_;
N=100;


[ InitialConditions.x0,InitialConditions.R0] = solveTime0Problem(Para,c,V,btild_1,s_);

parfor n=1:N
rHist0=rand(NumSim,1);
[SimData(n)]=tempRunSimulationsFromT1AltThetaShocksUsingPolicyRules(PolicyFunctions,InitialConditions,NumSim,Para,rHist0);
end
save('MCMC.mat','SimData')
clear SimData
clear InitialConditions



btild_1Grid=linspace(-2.5,2.5,NumInit);
s_=1;
rHist0=rand(NumSim,1);
tempInitialConditions.s0=0;
tempInitialConditions.x0=0;
tempInitialConditions.R0=0;

for n=1:NumInit
    
InitialConditions(n)=tempInitialConditions
end

parfor n=1:NumInit
    disp(n)
[ x0,R0] = solveTime0Problem(Para,c,V,btild_1Grid(n),s_);
InitialConditions(n).s0=s_;
InitialConditions(n).x0=x0;
InitialConditions(n).R0=R0;
[SimData(n)]=tempRunSimulationsFromT1AltThetaShocksUsingPolicyRules(PolicyFunctions,InitialConditions(n),NumSim,Para,rHist0);

end

save('DiffStarts.mat','SimData')

clear SimData
clear InitialConditions

NumSim=NumSim*10;
rHist0=rand(NumSim,1);
btild_1=-2;
s_=1;
InitialConditions.s0=s_;
[ InitialConditions.x0,InitialConditions.R0] = solveTime0Problem(Para,c,V,btild_1,s_);
SimData=tempRunSimulationsFromT1AltThetaShocksUsingPolicyRules(PolicyFunctions,InitialConditions,NumSim,Para,rHist0);
save('LongSample.mat','SimData')
break

load('MCMC2.mat','SimData')
% x
figure()
plot([SimData.xHist],':k')
hold on
plot(mean([SimData.xHist],2),'r','LineWidth',2)
%axis([1 length(SimData(1).xHist) Para.xMin Para.xMax])
% R
figure()
plot([SimData.RHist],':k')
hold on
plot(mean([SimData.RHist],2),'r','LineWidth',2)
%axis([1 length(SimData(1).RHist) Para.RMin Para.RMax])
a=2
figure()
Data=-[SimData.btildHist]./[SimData.YHist];
plot(Data(a:end,:),':k')
hold on
plot(mean([Data(a:end,:)],2),'r','LineWidth',2)
axis([1 length(Data(a:end,1)) 0 1])
print('-dpng','MCMCDebtGDP.png')

load('DiffStarts2.mat','SimData')
a=2
figure()
Data=-[SimData.btildHist]./[SimData.YHist];
plot(Data(a:end,:),'k')
print('-dpng','DiffStartDebtGDP.png')

load('LongSample.mat','SimData')
a=2
figure()
Data=-SimData.btildHist./SimData.YHist;
plot(Data(a:end,:),'k','LineWidth',2)
print('-dpng','LongSampleDebtGDP.png')

load('LongSampleZero.mat','SimData')
a=2
figure()
Data=-SimData.btildHist./SimData.YHist;
plot(Data(a:end,:),'k','LineWidth',2)
print('-dpng','LongSampleDebtGDPZero.png')
