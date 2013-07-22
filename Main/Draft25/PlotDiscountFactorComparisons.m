close all
opengl software
clear all
NumSim=25;
load rHist0.mat
rHist0=rHist0(1:NumSim)
%rHist0=rand(NumSim,1);
btild_1=-2;
s_=1;
InitialConditions.s0=s_;
%load ~/Golosov-Sargent/Data/temp/Simulations.mat
% IID - no beta shocks



load ~/Golosov-Sargent/Data/Draft/cTFPIneqLargerBetaShocks.mat
%[PolicyFunctions,errorPolicyFunctions]=GetPolicyRuleApproximations(Para,c,V,35,35,PolicyRulesStore,'spli',25,25,3,domain);
%save('~/Golosov-Sargent/Data/Draft/cTFPIneqBetaShocks.mat' , 'c','ErrorInSupNorm','cdiff','IndxSolved','IndxUnSolved','PolicyRulesStore','VNew','domain','Para','V','xhat','Coeff_xhat','Rhat','Coeff_Rhat','PolicyFunctions','errorPolicyFunctions');
[ InitialConditions.x0,InitialConditions.R0] = solveTime0Problem(Para,c,V,btild_1,s_);
[SimData(1)]=RunSimulationsFromT1AltThetaShocksUsingPolicyRules(PolicyFunctions,InitialConditions,NumSim,Para,rHist0);

load ~/Golosov-Sargent/Data/Draft/cTFPIneq.mat
[ InitialConditions.x0,InitialConditions.R0] = solveTime0Problem(Para,c,V,btild_1,s_);
[SimData(2)]=RunSimulationsFromT1AltThetaShocksUsingPolicyRules(PolicyFunctions,InitialConditions,NumSim,Para,rHist0);


figure()
subplot(2,2,1)
plot([SimData(1).TauHist(2:end)],'k','LineWidth',2.5)
hold on
plot([SimData(2).TauHist(2:end)],':k','LineWidth',2.5)
title('Labor Tax')
xlabel('time')

 bbT{1}=0:1;
Inx=find(SimData(1).sHist(2:end)<2);
    for i=2:length(Inx)-1
                bbT{i}=Inx(i)-.5:min(Inx(i)+.5,NumSim); 
    end
axis tight
ShadePlotForEmpahsis( bbT,'r',.05);  
subplot(2,2,2)
plot([SimData(1).TransHist(2:end)],'k','LineWidth',2.5)
hold on
plot([SimData(2).TransHist(2:end)],':k','LineWidth',2.5)
title('Transfers')
xlabel('time')

 bbT{1}=0:1;
Inx=find(SimData(1).sHist(2:end)<2);
    for i=2:length(Inx)
        bbT{i}=Inx(i)-.5:min(Inx(i)+.5,NumSim);
    end
axis tight
ShadePlotForEmpahsis( bbT,'b',.05);  



subplot(2,2,3)
plot(-[SimData(1).btildHist(2:end)],'k','LineWidth',2.5)
hold on
plot(-[SimData(2).btildHist(2:end)],':k','LineWidth',2.5)
 bbT{1}=0:1;

 Inx=find(SimData(1).sHist(2:end)<2);
    for i=2:length(Inx)-1
                bbT{i}=Inx(i)-.5:min(Inx(i)+.5,NumSim); 
    end
axis tight
ShadePlotForEmpahsis( bbT,'b',.05);  
title('Debt')
h = get(gca, 'title');

xlabel('time')

subplot(2,2,4)
plot([SimData(1).IntHist(2:end)-1],'k','LineWidth',2.5)
hold on
plot([SimData(2).IntHist(2:end)-1],':k','LineWidth',2.5)

 bbT{1}=0:1;
Inx=find(SimData(1).sHist(2:end)<2);
    for i=2:length(Inx)-1
                bbT{i}=Inx(i)-.5:min(Inx(i)+.5,NumSim); 
    end
axis tight
ShadePlotForEmpahsis( bbT,'b',.05);  
title('Interest Rates')
xlabel('time')

set(findall(gcf,'type','text'),'FontSize',14,'fontWeight','bold')