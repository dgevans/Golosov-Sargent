% This file plots the steady state decomposition of government financing

% load the case
clear all
close all
rootDirTemp=pwd;
rootDir=rootDirTemp(1:end-length(['/SteadyStateCode'])); % get root directory
run([rootDir '/Main/SetPath.m'])
InitData=load([rootDir '/Data/temp/csigmaLow.mat']);
texpath=[rootDir 'Tex/'];
Para=InitData.Para;
plotpath=[rootDir '/Graphs']
% Define the comparative static experiment
thetaRatioMin=2;
thetaRatioMax=5;
thetaRatioGridSize=25;
thetaRatioGrid=linspace(thetaRatioMin,thetaRatioMax,thetaRatioGridSize);
for i=1:thetaRatioGridSize
    Para=InitData.Para;
    Para.theta_2=1;
    Para.theta_1=Para.theta_2*thetaRatioGrid(i);
    Param(i)=Para;
end
ParamName='theta_1'
ParamBMValue=InitData.Para.theta_1
ParamGrid=thetaRatioGrid
GetComparitiveStaticsPlot(Param,ParamName,ParamGrid,ParamBMValue)



alpha1Min=.15;
alpha1Max=.85;
alpha1GridSize=25;
alpha1Grid=linspace(alpha1Min,alpha1Max,alpha1GridSize);
for i=1:alpha1GridSize
    Para=InitData.Para;
    Para.alpha_1=alpha1Grid(i);
Para.alpha_2=1-alpha1Grid(i);
    Param(i)=Para;


end

ParamName='alpha_1'
ParamBMValue=InitData.Para.alpha_1
ParamGrid=alpha1Grid
GetComparitiveStaticsPlot(Param,ParamName,ParamGrid,ParamBMValue)


sigmaMin=.8;
sigmaMax=3;
sigmaGridSize=25;
sigmaGrid=linspace(sigmaMin,sigmaMax,sigmaGridSize);

for i=1:sigmaGridSize
    Para=InitData.Para;
    Para.sigma=sigmaGrid(i);
    Param(i)=Para;
end


ParamName='sigma_1'
ParamBMValue=InitData.Para.sigma
ParamGrid=sigmaGrid
GetComparitiveStaticsPlot(Param,ParamName,ParamGrid,ParamBMValue,plotpath)
