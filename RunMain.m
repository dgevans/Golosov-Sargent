 clc
 clear all
 close all

SetParaStruc
Para.datapath=['Data/Productivity/'];
mkdir(Para.datapath)
casename='Productivity';
Para.StoreFileName=['c' casename '.mat'];
CoeffFileName=[Para.datapath Para.StoreFileName];
 
 %  --- SOLVE THE BELLMAN EQUATION --------------------------------------
 % test run 
Para.Niter=25;
RGrid.RMin=2.5;
RGrid.RMax=2.8;
LoadIndx=MainBellman(Para,RGrid);
NumIter=LoadIndx;
while NumIter < 20
    %InitData = load(CoeffFileName);
InitData=load([Para.datapath 'c_' num2str(LoadIndx) '.mat']);
RGrid.RMax=min(InitData.x_state(InitData.IndxUnSolved,2))*.95;
RGrid.RMin=2.5;
LoadIndx=MainBellman(Para,RGrid,InitData);
NumIter=NumIter+LoadIndx;
end

Para.Niter=5;
MainBellman(Para,RGrid);

%-- Simulate the MODEL -------------------------------------------------
NumSim=250;
sHist0=round(rand(NumSim,1))+1;


K=1;

ex(1).casename='Productivity'; % benchmark calibrations high alpha1



for ctrb=1:K
CoeffFileName=['Data/Productivity/c' ex(ctrb).casename '.mat'];
Sol=load(CoeffFileName);
Param(ctrb)=Sol.Para;
c10guess=1;
c20guess=1/Param(ctrb).RMax;
end

for ctrb=1:K
  CoeffFileName=['Data/Productivity/c' ex(ctrb).casename '.mat'];
c10guess=1;
c20guess=1/Param(ctrb).RMax;
[sHist,u2btildHist,RHist,TauHist,YHist,TransHist,...
          btildHist,c1Hist,c2Hist,l1Hist,l2Hist,IntHist,...
          IncomeFromAssets_Agent1Hist,AfterTaxWageIncome_Agent1Hist,...
          AfterTaxWageIncome_Agent2Hist,TransDiffHist,...
          LaborTaxAgent1DiffHist,LaborTaxAgent2DiffHist,DebtDiffHist,...
          GiniCoeffHist]=RunSimulations(CoeffFileName,NumSim,Para,sHist0);
end

save( [Para.datapath 'SimDataParallelProductivity.mat'],'sHist',...
       'u2btildHist','RHist','TauHist','YHist','TransHist',...
       'btildHist','c1Hist','c2Hist','l1Hist','l2Hist','Para','IntHist',...
       'AfterTaxWageIncome_Agent1Hist','AfterTaxWageIncome_Agent2Hist',...
       'IncomeFromAssets_Agent1Hist','TransDiffHist',...
       'LaborTaxAgent1DiffHist','LaborTaxAgent2DiffHist','DebtDiffHist',...
       'GiniCoeffHist')
    
   
   
   
%  % -- PLOT DIAGNOSTICS -----------------------------------------
close all
clear all
clc
SimTitle={'Productivity'};
SimDataPath= 'Data/Productivity/SimDataParallelProductivity.mat';
SimPlotPath='Graphs/Productivity/';
mkdir(SimPlotPath)
SimTexPath='Tex/Productivity/';
mkdir(SimTexPath)
PlotParallelSimulationsCommonShocks(SimDataPath,SimTexPath,SimPlotPath,SimTitle)


 Para.datapath=['Data/'];
 Para.StoreFileName=['c_4.mat'];
 
 GetPlotsForFinalSolution(Para)
