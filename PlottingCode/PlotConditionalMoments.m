% This file plots the transition dynamics in the XR space
clear all
clc
close all
rootDirTemp=pwd;
rootDir=rootDirTemp(1:end-length(['/PlottingCode'])); % get root directory
run([rootDir '/Main/SetPath.m'])

%ex(1).casename='sigmaLow';
ex(1).casename='inequality';
i=1
% load the coeffecients
BellmanData=load([rootDir sl 'Data' sl 'temp' sl 'c' ex(i).casename '.mat']);
 BellmanData.Para.plotpath=[rootDir sl 'Graphs/'];
% load simulation data
 BellmanData.Para.datapath=[rootDir sl 'Data/temp/'];
 BellmanData.Para.StoreFileName=['c' ex(i).casename '.mat'];
 BellmanData.Para.flagPlot2PeriodDrifts=0;
 ComputeConditionalMoments(BellmanData.Para)
 close all
 GetPlotsForFinalSolution(BellmanData.Para)
 close all
clc
SimTitle{1}=' inequality';
SimDataPath= 'C:\Users\Anmol\Dropbox\2011RA\FiscalPolicy\OrganizedCode\Golosov-Sargent\Data\temp\SimDataThetaShocks.mat';
SimPlotPath=BellmanData.Para.plotpath;
mkdir(SimPlotPath)
SimTexPath='Tex/';
mkdir(SimTexPath)
PlotParallelSimulationsCommonShocks(SimDataPath,SimTexPath,SimPlotPath,SimTitle,BellmanData.Para)
