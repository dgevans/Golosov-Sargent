% This file plots the transition dynamics in the XR space
clear all
clc
close all

load ('~/Golosov-Sargent/Data/temp/cTFPHighInequalityBeta.mat');
%SimTitle{1}=' : expenditure shocks';
SimTitle{1}=' : prod+ineq+beta shocks';
SimDataPath= '~/Golosov-Sargent/Data/temp/SimulationTFPHighIneqBeta.mat'
SimPlotPath='/home/anmol/Dropbox/2011RA/FiscalPolicy/GolosovProjectCode/SteadyStateNotes/Graphs/ThetaShocks/';
mkdir(SimPlotPath)
%SimTexPath='~/Golosov-Sargent/\Tex\GShocks\';
SimTexPath='/home/anmol/Dropbox/2011RA/FiscalPolicy/GolosovProjectCode/SteadyStateNotes/Tex/ThetaShocks/';
mkdir(SimTexPath)
Para.BigT=3000;
PlotParallelSimulationsCommonShocks(SimDataPath,SimTexPath,SimPlotPath,SimTitle,Para)
