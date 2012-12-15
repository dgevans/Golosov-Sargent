plotpath='C:\Users\Anmol\Dropbox\2011RA\FiscalPolicy\GolosovProjectCode\Tom Example\BGP\Graphs\War\';

mkdir(plotpath)
load data\Calibration\SimDataParallelWar.mat
% figure()
plot(u2btildHist)
xlabel('x_t')
ylabel('t')
title('x')
print(gcf,'-dpng',[plotpath 'LongSimulationX.png'])
figure()
plot(RHist)
xlabel('R_t')
ylabel('t')
title('R')
print(gcf,'-dpng',[plotpath 'LongSimulationR.png'])


%  % -- PLOT simulations -----------------------------------------
close all
clear all
clc
SimTitle{1}='War';
SimDataPath= 'Data/Calibration/SimDataParallelWar.mat';
SimPlotPath='C:\Users\Anmol\Dropbox\2011RA\FiscalPolicy\GolosovProjectCode\Tom Example\BGP\Graphs\War\';
mkdir(SimPlotPath)
SimTexPath='C:\Users\Anmol\Dropbox\2011RA\FiscalPolicy\GolosovProjectCode\Tom Example\BGP\Tex\War\';
mkdir(SimTexPath)
PlotParallelSimulationsCommonShocks(SimDataPath,SimTexPath,SimPlotPath,SimTitle)