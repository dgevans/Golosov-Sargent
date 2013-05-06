% This file plots the transition dynamics in the XR space
clear all
clc
close all

load ('~/Golosov-Sargent/Data/temp/cHighIneq.mat');
%SimTitle{1}=' : expenditure shocks';
SimTitle{1}=' : prod+ineq shocks';
SimDataPath= '~/Golosov-Sargent/Data/temp/SimulationHighIneq'
SimPlotPath='~/Golosov-Sargent/Plots/PrShocks/';
mkdir(SimPlotPath)
%SimTexPath='~/Golosov-Sargent/\Tex\GShocks\';
SimTexPath='~/Golosov-Sargent/TexPrShocks/';
mkdir(SimTexPath)
Para.BigT=25000;
PlotParallelSimulationsCommonShocks(SimDataPath,SimTexPath,SimPlotPath,SimTitle,Para)
