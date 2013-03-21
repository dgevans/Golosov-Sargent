
plotpath='C:\Users\Anmol\Dropbox\2011RA\FiscalPolicy\GolosovProjectCode\Tom Example\BGP\Graphs\AltStarts\';

mkdir(plotpath)
load data\Calibration\SimDataStationaryStart.mat
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
SimTitle{1}='AltStarts';
SimDataPath= 'Data/Calibration/SimDataStationaryStart.mat';
SimPlotPath='C:\Users\Anmol\Dropbox\2011RA\FiscalPolicy\GolosovProjectCode\Tom Example\BGP\Graphs\AltStarts\';
mkdir(SimPlotPath)
SimTexPath='C:\Users\Anmol\Dropbox\2011RA\FiscalPolicy\GolosovProjectCode\Tom Example\BGP\Tex\AltStarts\';
mkdir(SimTexPath)
PlotParallelSimulationsCommonShocks(SimDataPath,SimTexPath,SimPlotPath,SimTitle)