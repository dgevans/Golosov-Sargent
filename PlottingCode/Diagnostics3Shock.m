clear all
clc
close all
rootDirTemp=pwd;
rootDir=rootDirTemp(1:end-length(['/PlottingCode'])); % get root directory
run([rootDir '/Main/SetPath.m'])

ex(1).casename='sigmalow'; 
ex(2).casename='sigmaMed';
ex(3).casename='sigmaHigh'

i=2
% load the coeffecients
BellmanData=load(['Data/temp/c' ex(i).casename '.mat']);
 BellmanData.Para.plotpath=[rootDir sl 'Graphs/3Shocks'];
 BellmanData.Para.datapath=[rootDir sl 'Data/temp/'];

 PlotValueFunction(BellmanData.Para)
 load([rootDir '/Data/temp/SimData3ShockSigmaMed.mat'])
 figure()
 subplot(1,2,1)
 plot(SimData.xHist,'k','LineWidth',2)
 xlabel('time')
 ylabel('x')
 title('Simulations for $\sigma =2$','Interpreter','Latex')
 subplot(1,2,2)
 plot(SimData.RHist,'k','LineWidth',2)
 xlabel('time')
 ylabel('\rho')
 title('Simulations for $\sigma =2$','Interpreter','Latex')