 clc
 clear all
 close all

 %matlabpool open
run('~/Golosov-Sargent/Main/SetPath')
load ('~/Golosov-Sargent/Data/temp/cTFPHighInequalityBeta.mat');
g=Para.g;
n1=Para.n1;
n2=Para.n2;
alpha_1=Para.alpha_1;
alpha_2=Para.alpha_2;
theta_1=Para.theta_1;
theta_2=Para.theta_2;
psi=Para.psi;
beta=Para.beta;
sigma=Para.sigma;
Para.U=@(c,l) UMix(c,l,Para);
NumSim=50000;
K=1;
ex(1).casename='b_{-1}=-1'; 
Para.saveSimPath= ['~/Golosov-Sargent/Data/temp/SimulationTFPHighIneqBeta.mat'];
rHist0 = rand(NumSim,1);
CoeffFileName='~/Golosov-Sargent/Data/temp/cTFPHighInequalityBeta.mat';
tempbeta=Para.beta;
Para.beta=mean(Para.beta);
[x, R,PolicyRule]=findSteadyState(0,2,Para);
Para.beta=tempbeta;
[x, R,PolicyRule]=findSteadyState(x,R,Para,PolicyRule);
% EXP
SimDataHighIneq=RunSimulationsFromT1AltThetaShocks(CoeffFileName,-1,R,NumSim,Para,rHist0);
