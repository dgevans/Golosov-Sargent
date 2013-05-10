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
ex(1).casename='HighAlpha1';

i=1
% LOAD COEFF
 BellmanData=load([rootDir sl 'Data' sl 'temp' sl 'c' ex(i).casename '.mat']);
 BellmanData.Para.plotpath=[rootDir sl 'Graphs/'];
 BellmanData.Para.datapath=[rootDir sl 'Data/temp/'];
 BellmanData.Para.StoreFileName=['c' ex(i).casename '.mat'];
 
 
 
 %% PLOT CONDITIONAL MOMENTS
 %BellmanData.Para.flagPlot2PeriodDrifts=0;
 %BellmanData.Para.flagComputeAutoCorr=1;
 %ComputeConditionalMoments(BellmanData.Para)

 
 
 % PLOT POLICYRULES
 
 
 close all
 %GetPlotsForFinalSolution(BellmanData.Para)
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


% SOLVE THE T-0 PROBLEM given btild(-1)
btild_1=-1;
s_=1;
c10guess=1;
c20guess=.5;
%disp('Computed V...Now solving V0(btild_1) where btild_1 is')
%disp(btild_1)
% c1 and c2 solve
options=optimset('Display','off');
[x,~,exitflagv0,~,~] = fminunc(@(x)  getValue0(x, btild_1,s_,BellmanData.Para,BellmanData.c,BellmanData.V),[ c10guess c20guess],options);

c10 = x(1);
c20 = x(2);
R0=(c10/c20)^(sigma);
TotalResources=(c10*n1+c20*n2+g(s_));
DenL2=theta_2*R0*n1+theta_2*n2;
l20=(TotalResources-theta_1*n1+ theta_2*n1*R0)/(DenL2);
l10= 1-(1-l20)*theta_2/theta_1*R0;
x0=-(c20-c10)*(psi*c20^(-sigma))-((l10/(1-l10))*R0-l20/(1-l20))*(1-psi)+btild_1*psi*c20^(-sigma);
R0=c20^(-sigma)/c10^(-sigma);

%NormDecomp=ComputeConditionalChangeBudgetConstraint(BellmanData.Para,BellmanData.c,BellmanData.V,1,BellmanData.domain,BellmanData.PolicyRulesStore,x0,R0);

%NormDecomp



NormDecomp=ComputeConditionalChangeBudgetConstraint(BellmanData.Para,BellmanData.c,BellmanData.V,1,BellmanData.domain,BellmanData.PolicyRulesStore,0,BellmanData.Para.RSS);

NormDecomp