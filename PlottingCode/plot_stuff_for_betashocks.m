% This file plots the transition dynamics in the XR space
clear all
%clc
close all

%% RUN THIS SCRIPT FROM THE PLOTTING CODE DIRECTORY
rootDir='~/Golosov-Sargent'
run([rootDir '/Main/SetPath.m'])

%ex(1).casename='sigmaLow';
%ex(1).casename='HighIneq';
%ex(1).casename='productivity';
ex(1).casename='symIneq';

i=1
% LOAD COEFF
 BellmanData=load([rootDir sl 'Data' sl 'temp' sl 'c' ex(i).casename '.mat']);
 BellmanData.Para.plotpath=[rootDir sl 'Graphs/'];
 BellmanData.Para.datapath=[rootDir sl 'Data/temp/'];
 BellmanData.Para.StoreFileName=['c' ex(i).casename '.mat'];
 
 
 
 %% PLOT CONDITIONAL MOMENTS
 BellmanData.Para.flagPlot2PeriodDrifts=0;
 BellmanData.Para.flagComputeAutoCorr=1;
 %ComputeConditionalMoments(BellmanData.Para)

 
 
 % PLOT POLICYRULES
 
 
 close all
%BellmanData.Para.U=@(c,l) UMix(c,l,BellmanData.Para)
%[ BellmanData.Para.xSS,BellmanData.Para.RSS,~ ] = findSteadyState( 0,mean(BellmanData.Para.RGrid),BellmanData.Para)
BellmanData.Para.xSS=0;
BellmanData.Para.RSS=mean(BellmanData.Para.RGrid);

 GetPlotsForFinalSolution(BellmanData.Para)
 close all
 % PLOT SIMULATIONS
clc
SimTitle{1}=' : inequality shocks';
SimDataPath= 'C:\Users\Anmol\Dropbox\2011RA\FiscalPolicy\OrganizedCode\Golosov-Sargent\Data\temp\SimDataThetaShocks.mat';
%SimTitle{1}=' : expenditure shocks';
%SimDataPath= 'C:\Users\Anmol\Dropbox\2011RA\FiscalPolicy\OrganizedCode\Golosov-Sargent\Data\temp\SimDataGShocksLongSample.mat';
%SimTitle{1}=' : productivity shocks';
%SimDataPath= 'C:\Users\Anmol\Dropbox\2011RA\FiscalPolicy\OrganizedCode\Golosov-Sargent\Data\temp\SimDataGShocksLongSampleProductivity.mat';


SimPlotPath=BellmanData.Para.plotpath;
mkdir(SimPlotPath)
SimTexPath='Tex/';
mkdir(SimTexPath)
BellmanData.Para.BigT=25000;
%PlotParallelSimulationsCommonShocks(SimDataPath,SimTexPath,SimPlotPath,SimTitle,BellmanData.Para)



g=BellmanData.Para.g;
n1=BellmanData.Para.n1;
n2=BellmanData.Para.n2;
alpha_1=BellmanData.Para.alpha_1;
alpha_2=BellmanData.Para.alpha_2;
theta_1=BellmanData.Para.theta_1;
theta_2=BellmanData.Para.theta_2;
psi=BellmanData.Para.psi;
beta=BellmanData.Para.beta;
sigma=BellmanData.Para.sigma;


%NormDecomp=ComputeConditionalChangeBudgetConstraint(BellmanData.Para,BellmanData.c,BellmanData.V,1,BellmanData.domain,BellmanData.PolicyRulesStore,x0,R0);

%NormDecomp

[x,R]=findSteadyState(0,3,BellmanData.Para)
NormDecomp=ComputeConditionalChangeBudgetConstraint(BellmanData.Para,BellmanData.c,BellmanData.V,1,BellmanData.domain,BellmanData.PolicyRulesStore,0,R);

NormDecomp