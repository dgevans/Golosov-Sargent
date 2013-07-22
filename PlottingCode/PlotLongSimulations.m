% This file plots the transition dynamics in the XR space
clear all
clc
close all
load('~/Golosov-Sargent/Data/temp/Simulations.mat');
load ('~/Golosov-Sargent/Data/temp/cTFPIneq.mat');
%SimTitle{1}=' : expenditure shocks';
SimTitle{1}=' : prod+ineq';
SimPlotPath='/home/anmol/Dropbox/2011RA/FiscalPolicy/GolosovProjectCode/Calibration/Graphs/NewCalibration/ThetaShocks/';
mkdir(SimPlotPath)
%SimTexPath='~/Golosov-Sargent/\Tex\GShocks\';
SimTexPath='/home/anmol/Dropbox/2011RA/FiscalPolicy/GolosovProjectCode/Calibration/Tex/NewCalibration/ThetaShocks/';
mkdir(SimTexPath)
Para.BigT=15000;
PlotParallelSimulationsCommonShocks(SimData(2),SimTexPath,SimPlotPath,SimTitle,Para)
figure()
subplot(2,1,1)
plot([SimData(1:4).xHist],'LineWidth',2)
legend('Constant beta','Low beta','Med beta','High beta')
ylabel('x_t ')
title('Marginal Utility Adjusted Difference Debt')

subplot(2,1,2)
plot([SimData(1:4).RHist],'LineWidth',2)
legend('Constant beta','Low beta','Med beta','High beta')
ylabel('\rho_t ')
title('Ratio of Marginal Utility')
