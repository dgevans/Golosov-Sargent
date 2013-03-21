plotpath='C:\Users\Anmol\Dropbox\2011RA\FiscalPolicy\GolosovProjectCode\Tom Example\BGP\Graphs\PExperiments\';

mkdir(plotpath)
load data\Calibration\SimDataParallelPCommonShocks.mat
figure()
plot(u2btildHist)
xlabel('x_t')
ylabel('t')
title('x')
legend('Ph/Pl<1','Ph/Pl=1','Ph/Pl>1')
print(gcf,'-dpng',[plotpath 'LongSimulationX.png'])
figure()
plot(RHist)
xlabel('R_t')
ylabel('t')
title('R')
legend('Ph/Pl<1','Ph/Pl=1','Ph/Pl>1')
print(gcf,'-dpng',[plotpath 'LongSimulationR.png'])


%  % -- PLOT simulations -----------------------------------------
close all
clear all
clc
SimTitle{1}='$Ph<Pl$';
SimTitle{2}='$Ph=Pl$';
SimTitle{3}='$Ph>Pl$';
SimDataPath= 'Data/Calibration/SimDataParallelPCommonShocks.mat';
SimPlotPath='C:\Users\Anmol\Dropbox\2011RA\FiscalPolicy\GolosovProjectCode\Tom Example\BGP\Graphs\PExperiments\';
mkdir(SimPlotPath)
SimTexPath='C:\Users\Anmol\Dropbox\2011RA\FiscalPolicy\GolosovProjectCode\Tom Example\BGP\Tex\PExperiments\';
mkdir(SimTexPath)
PlotParallelSimulationsCommonShocks(SimDataPath,SimTexPath,SimPlotPath,SimTitle)

% Plot Drifts and Volatilities
plotpath='C:\Users\Anmol\Dropbox\2011RA\FiscalPolicy\GolosovProjectCode\Tom Example\BGP\Graphs\PExperiments\';
K=3;
ex(1).casename='PhLow'; 
ex(2).casename='PhMed'; 
ex(3).casename='PhHigh';
for i=2:2
Domain.xBounds=[1.4 1.55];
Domain.RBounds=[2.75 2.9];
 Para.datapath=['Data/Calibration/'];
 Para.plotpath=[plotpath ex(i).casename '/'];
 Para.StoreFileName=['c' ex(i).casename '.mat'];
 GetPlotsForFinalSolution(Para,Domain)
end

