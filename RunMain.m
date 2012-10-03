 clc
 clear all
 close all

SetParaStruc
Para.datapath=['Data/temp/'];
mkdir(Para.datapath)
casename='Productivity';
Para.StoreFileName=['c' casename '.mat'];
CoeffFileName=[Para.datapath Para.StoreFileName];
 
 %  --- SOLVE THE BELLMAN EQUATION --------------------------------------
 % test run 
Para.Niter=150;
RGrid.RMin=2.5;
RGrid.RMax=3.2;
Para.datapath=['Data/temp/Productivity/'];
% LoadIndx=MainBellman(Para,RGrid);
% NumIter=LoadIndx;
% while NumIter < 200
%     %InitData = load(CoeffFileName);
% InitData=load([Para.datapath 'c_' num2str(LoadIndx) '.mat']);
% RGrid.RMax=min(InitData.x_state(InitData.IndxUnSolved,2))*.95;
% RGrid.RMin=2.5;
% LoadIndx=MainBellman(Para,RGrid,InitData);
% NumIter=NumIter+LoadIndx;
% end

% --- SOLVE THE BEllMAN FOR INEQUALITY SHOCKS -----------
meantheta = mean([theta_1_low,theta_2_low]);
Para.theta_1 = [theta_1_low;theta_1_low-0.05*meantheta];
Para.theta_2 = [theta_2_low;theta_2_low+0.05*meantheta];
Para.datapath=['Data/temp/Inequality/'];
mkdir(Para.datapath)
casename='Inequality';
Para.StoreFileName=['c' casename '.mat'];
CoeffFileName=[Para.datapath Para.StoreFileName];

LoadIndx=MainBellman(Para,RGrid);
NumIter=LoadIndx;
while NumIter < 100
    %InitData = load(CoeffFileName);
InitData=load([Para.datapath 'c_' num2str(LoadIndx) '.mat']);
RGrid.RMax=min(InitData.x_state(InitData.IndxUnSolved,2))*.95;
RGrid.RMin=2.5;
LoadIndx=MainBellman(Para,RGrid,InitData);
NumIter=NumIter+LoadIndx;
end

%-- Simulate the MODEL -------------------------------------------------
NumSim=10000;
sHist0=randi(2,NumSim,1);


K=2;

ex(1).casename='Productivity'; % benchmark calibrations high alpha1
ex(2).casename='Inequality';


for ctrb=1:K
CoeffFileName=['Data/temp/' ex(ctrb).casename '/c' ex(ctrb).casename '.mat'];
Sol=load(CoeffFileName);
Param(ctrb)=Sol.Para;
c10guess=1;
c20guess=1/Param(ctrb).RMax;
end

parfor ctrb=1:K
  CoeffFileName=['Data/temp/' ex(ctrb).casename '/c' ex(ctrb).casename '.mat'];
  
  [sHist(:,ctrb),u2btildHist(:,ctrb),RHist(:,ctrb),TauHist(:,ctrb),YHist(:,ctrb),TransHist(:,ctrb),...
          btildHist(:,ctrb),c1Hist(:,ctrb),c2Hist(:,ctrb),l1Hist(:,ctrb),l2Hist(:,ctrb),IntHist(:,ctrb),...
          IncomeFromAssets_Agent1Hist(:,ctrb),AfterTaxWageIncome_Agent1Hist(:,ctrb),...
          AfterTaxWageIncome_Agent2Hist(:,ctrb),TransDiffHist(:,ctrb),...
          LaborTaxAgent1DiffHist(:,ctrb),LaborTaxAgent2DiffHist(:,ctrb),DebtDiffHist(:,ctrb),...
          GiniCoeffHist(:,ctrb),theta_1Hist(:,ctrb),theta_2Hist(:,ctrb)]=RunSimulations(CoeffFileName,NumSim,Para,sHist0);
end
Para.datapath=['Data/temp/'];

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
SimTitle={'Productivity','Inequality'};
SimDataPath= 'Data/temp/SimDataParallel.mat';
SimPlotPath='Graphs/Parallel/';
mkdir(SimPlotPath)
SimTexPath='Tex/Parallel/';
mkdir(SimTexPath)
PlotParallelSimulationsCommonShocks(SimDataPath,SimTexPath,SimPlotPath,SimTitle)


 Para.datapath=['Data/temp/Inequality/'];
 Para.StoreFileName=['cInequality.mat'];
 
 GetPlotsForFinalSolution(Para)
