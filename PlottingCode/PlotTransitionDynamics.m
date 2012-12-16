% This file plots the transition dynamics in the XR space
clear all
clc
close all
rootDirTemp=pwd;
rootDir=rootDirTemp(1:end-length(['/PlottingCode'])); % get root directory
run([rootDir '/Main/SetPath.m'])

ex(1).casename='sigmaLow'; 
ex(2).casename='sigmaMed';
ex(3).casename='sigmaHigh'

i=1
% load the coeffecients
BellmanData=load(['Data/temp/c' ex(i).casename '.mat']);
 BellmanData.Para.plotpath=[rootDir sl 'Graphs/'];
% load simulation data
load (['Data/temp/SimDataParallelCommonShocks.mat'])
figure()
subplot(2,1,1)
plot(xHist(:,i),'k','LineWidth',2)
xlabel('t')
ylabel('x_t')
axis([ 1 length(xHist(:,i))  min(xHist(:,i))  max(xHist(:,i))*1.2])
subplot(2,1,2)
plot(RHist(:,1),'k','LineWidth',2)
xlabel('t')
ylabel('$\rho_t$','Interpreter','Latex')
print('-dpng',[rootDir sl 'Graphs\LongSimulationsXR.png'])
% Checking against david's code
[ xSS,RSS,PolicyRule ] = findSteadyState( 0,mean(BellmanData.Para.RGrid),BellmanData.Para);
SimulationData.xHist=xHist(:,i);
SimulationData.RHist=RHist(:,i);
SimulationData.SampleFrequency=500;
SimulationData.InitialState_x=xHist(1,i);
SimulationData.InitialState_R=RHist(1,i);
SimulationData.SteadyState_x=xSS;
SimulationData.SteadyState_R=RSS;
GetTransitionDynamics(BellmanData,SimulationData)
GetPlotsForFinalSolution3D(BellmanData)
 BellmanData.Para.datapath=[rootDir sl 'Data/temp/'];
 BellmanData.Para.StoreFileName=['c' ex(i).casename '.mat'];
 BellmanData.Para.flagPlot2PeriodDrifts=0;
 GetPlotsForFinalSolution(BellmanData.Para)
 C={'k','r','g'}
 figure()
 for i=1:3
     BellmanData=load(['BellmanData.Para.datapath' 'c' ex(i).casename '.mat']);

[ xSS(i),RSS(i),PolicyRule ] = findSteadyState( 0,mean(BellmanData.Para.RGrid),BellmanData.Para);

 plot(xHist(:,i),RHist(:,i),C{i},'LineWidth',2)
 hold on
 
 end
 
 legend('\sigma=1','\sigma=2','\sigma=3')
 
 hold on
 for i=1:3
 scatter(xSS,RSS,C{i})
 hold on
 end
 xlabel('x')
 ylabel('$\rho$','Interpreter','Latex')
 