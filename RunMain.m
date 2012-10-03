 clc
 clear all
 close all

SetParaStruc
Para.datapath=['Data/temp/Productivity/'];
mkdir(Para.datapath)
casename='Productivity';
Para.StoreFileName=['c' casename '.mat'];
CoeffFileName=[Para.datapath Para.StoreFileName];
 
 %  --- SOLVE THE BELLMAN EQUATION --------------------------------------
 % test run 
Para.Niter=250;
RGrid.RMin=2.5;
RGrid.RMax=3.5;
LoadIndx=MainBellman(Para,RGrid);
NumIter=LoadIndx;
while NumIter < 200
    %InitData = load(CoeffFileName);
InitData=load([Para.datapath 'c_' num2str(LoadIndx) '.mat']);
RGrid.RMax=min(InitData.x_state(InitData.IndxUnSolved,2))*.95;
RGrid.RMin=2.5;
LoadIndx=MainBellman(Para,RGrid,InitData);
NumIter=NumIter+LoadIndx;
end

%-- Simulate the MODEL -------------------------------------------------
NumSim=10000;
sHist0=round(rand(NumSim,1))+1;


K=1;

ex(1).casename='Productivity'; % benchmark calibrations high alpha1



for ctrb=1:K
CoeffFileName=['Data/temp/Productivity/c' ex(ctrb).casename '.mat'];
Sol=load(CoeffFileName);
Param(ctrb)=Sol.Para;
c10guess=1;
c20guess=1/Param(ctrb).RMax;
end

for ctrb=1:K
  CoeffFileName=['Data/temp/Productivity/c' ex(ctrb).casename '.mat'];
  
  [sHist(:,ctrb),u2btildHist(:,ctrb),RHist(:,ctrb),TauHist(:,ctrb),YHist(:,ctrb),TransHist(:,ctrb),...
          btildHist(:,ctrb),c1Hist(:,ctrb),c2Hist(:,ctrb),l1Hist(:,ctrb),l2Hist(:,ctrb),IntHist(:,ctrb),...
          IncomeFromAssets_Agent1Hist(:,ctrb),AfterTaxWageIncome_Agent1Hist(:,ctrb),...
          AfterTaxWageIncome_Agent2Hist(:,ctrb),TransDiffHist(:,ctrb),...
          LaborTaxAgent1DiffHist(:,ctrb),LaborTaxAgent2DiffHist(:,ctrb),DebtDiffHist(:,ctrb),...
          GiniCoeffHist(:,ctrb),theta_1Hist(:,ctrb),theta_2Hist(:,ctrb)]=RunSimulations(CoeffFileName,NumSim,Para,sHist0);
end

save( [Para.datapath 'SimDataParallelProductivity.mat'],'sHist',...
       'u2btildHist','RHist','TauHist','YHist','TransHist',...
       'btildHist','c1Hist','c2Hist','l1Hist','l2Hist','Para','IntHist',...
       'AfterTaxWageIncome_Agent1Hist','AfterTaxWageIncome_Agent2Hist',...
       'IncomeFromAssets_Agent1Hist','TransDiffHist',...
       'LaborTaxAgent1DiffHist','LaborTaxAgent2DiffHist','DebtDiffHist',...
       'GiniCoeffHist','theta_1Hist','theta_2Hist')
    
   
   
   
%  % -- PLOT DIAGNOSTICS -----------------------------------------
close all
clear all
clc
SimTitle={'Productivity'};
SimDataPath= 'Data/temp/Productivity/SimDataParallelProductivity.mat';
SimPlotPath='Graphs/Productivity/';
mkdir(SimPlotPath)
SimTexPath='Tex/Productivity/';
mkdir(SimTexPath)
% PlotParallelSimulationsCommonShocks(SimDataPath,SimTexPath,SimPlotPath,SimTitle)
% 
% 
%  Para.datapath=['Data/temp/Productivity/'];
%  Para.StoreFileName=['c_25.mat'];
%  
%  GetPlotsForFinalSolution(Para)
