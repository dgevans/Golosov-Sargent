% This file plots the transition dynamics in the XR space
clear all
clc
close all

%% RUN THIS SCRIPT FROM THE PLOTTING CODE DIRECTORY
rootDirTemp=pwd;
rootDir=rootDirTemp(1:end-length(['/PlottingCode'])); % get root directory
run([rootDir '/Main/SetPath.m'])

ex(1).casename='sigmaLow';
%ex(1).casename='LowParetoAgent1';
%ex(1).casename='inequality';
i=1
% LOAD COEFF
 BellmanData=load([rootDir sl 'Data' sl 'temp' sl 'c' ex(i).casename '.mat']);
 BellmanData.Para.plotpath=[rootDir sl 'Graphs/'];
 BellmanData.Para.datapath=[rootDir sl 'Data/temp/'];
 BellmanData.Para.StoreFileName=['c' ex(i).casename '.mat'];
 
 
 
 %% PLOT CONDITIONAL MOMENTS
 BellmanData.Para.flagPlot2PeriodDrifts=0;
 BellmanData.Para.flagComputeAutoCorr=1;
 ComputeConditionalMoments(BellmanData.Para)

 
 
 % PLOT POLICYRULES
 
 
 close all
 GetPlotsForFinalSolution(BellmanData.Para)
 close all
 
 % PLOT SIMULATIONS
clc
%SimTitle{1}=' : inequality shocks';
%SimDataPath= 'C:\Users\Anmol\Dropbox\2011RA\FiscalPolicy\OrganizedCode\Golosov-Sargent\Data\temp\SimDataThetaShocks.mat';
SimTitle{1}=' : expenditure shocks';
SimDataPath= 'C:\Users\Anmol\Dropbox\2011RA\FiscalPolicy\OrganizedCode\Golosov-Sargent\Data\temp\SimDataGShocksDiffStarts.mat';
SimPlotPath=BellmanData.Para.plotpath;
mkdir(SimPlotPath)
SimTexPath='Tex/';
mkdir(SimTexPath)
PlotParallelSimulationsCommonShocks(SimDataPath,SimTexPath,SimPlotPath,SimTitle,BellmanData.Para)
