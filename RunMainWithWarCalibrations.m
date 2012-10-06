 clc
 clear all
 close all

 %{ 
This file solves the G-S economy with BGP preferences of the form
 psi.log[c]+(1-psi)log[1-l] with following calibrations

1.] The ratio of productivities is 3.3 and the low productivity is
 normalized to 1
2.] psi is choosen to get a average FE of about 0.5
3.] pareto wts are such that the no-shock problem gives a tax rate of about
    20 percent
4.] Government expenditures are about 11 and 13 percent of the output
5.] beta =0.9
%}
 
 
% Set up the para structure  
SetParaStruc
theta_1=3.3; % theta high
theta_2=1;  % theta low
g_l_y=.11; % g low
g_h_y=.13; % g high
n1=1;
n2=1;
tau=.2;
g_Y=mean([g_l_y g_h_y]);
AvfFETarget=.5;
x=fsolve(@(x) GetCalibrationFrischElasticity (x,AvfFETarget,theta_1,theta_2,tau,g_Y,n1,n2), [1 1 ]);
gamma=x(1);
Y=x(2);
g=g_Y*Y;
psi=1/(1+gamma);
beta=.9;
alpha_1=0.75;
alpha_2=1-alpha_1;
Para.n1=n1;
Para.n2=n2;
alpha_1=alpha_1*Para.n1;
alpha_2=alpha_2*Para.n2;
Para.beta=.9;
Para.alpha_1=alpha_1;
Para.alpha_2=alpha_2;
Para.psi=psi;
Para.g=[g_l_y g_h_y]*Y;

% WAR CALIBRATION
gMean=sum(Para.P(1,:).*Para.g);
Para.g(2)=2*Para.g(1);
NewPh=(gMean-Para.g(1))./(Para.g(2)-Para.g(1));
Para.P=[1-NewPh NewPh;1-NewPh NewPh];



Para.theta_1=theta_1;
Para.theta_2=theta_2;
Para.btild_1=0;
Para.alpha_1=alpha_1;
Para.alpha_2=alpha_2;
Para.datapath=['Data/temp/War/'];
mkdir(Para.datapath)
casename='War';
Para.StoreFileName=['c' casename '.mat'];
CoeffFileName=[Para.datapath Para.StoreFileName];
 
% % WAR CALIBRATION

HLRGridSize=1;
HLRMin=2;
HLRMax=2;
HLRGrid=linspace(HLRMin,HLRMax,HLRGridSize);

HighLowRatio=HLRGrid(1);
gMean=sum(Para.P(1,:).*Para.g);
Para.g(2)=HighLowRatio*Para.g(1);
NewPh=(gMean-Para.g(1))./(Para.g(2)-Para.g(1));
Para.P=[1-NewPh NewPh;1-NewPh NewPh];
Para.Niter=150;
RGrid.RMin=2.5;
RGrid.RMax=4;
LoadIndx=MainBellman(Para,RGrid);
InitData=load([Para.datapath 'c_' num2str(LoadIndx) '.mat']);
% 
% % g_h>g_l
% for i =2:HLRGridSize
%     Para.datapath=['Data/temp/War/'];
% mkdir(Para.datapath)
% casename='War';
% Para.StoreFileName=['c' casename num2str(i) '.mat'];
% CoeffFileName=[Para.datapath Para.StoreFileName];
% 
% HighLowRatio=HLRGrid(i);
% gMean=sum(Para.P(1,:).*Para.g);
% Para.g(2)=HighLowRatio*Para.g(1);
% NewPh=(gMean-Para.g(1))./(Para.g(2)-Para.g(1));
% Para.P=[1-NewPh NewPh;1-NewPh NewPh];
% Para.Niter=50;
% RGrid.RMin=2.5;
% RGrid.RMax=4;
% LoadIndx=MainBellman(Para,RGrid);
% NumIter=LoadIndx;
% while (NumIter < Para.Niter*.9 && RGrid.RMax>RGrid.RMin)
% InitData=load([Para.datapath 'c_' num2str(LoadIndx) '.mat']);
% RGrid.RMax=min(InitData.x_state(InitData.IndxUnSolved,2))*.98;
% RGrid.RMin=2.5;
% LoadIndx=MainBellman(Para,RGrid,InitData);
% NumIter=NumIter+LoadIndx;
% end
% if LoadIndx<Para.Niter*.8
%     break;
% end
% end
% 
% Para.Niter=150;
% MainBellman(Para,RGrid,InitData);

%-- Simulate the MODEL -------------------------------------------------
NumSim=10000;
sHist0=round(rand(NumSim,1))+1;


K=1;

ex(1).casename='War'; % benchmark calibrations high alpha1



for ctrb=1:K
CoeffFileName=['Data/temp/War/c' ex(ctrb).casename '.mat'];
Sol=load(CoeffFileName);
Param(ctrb)=Sol.Para;
c10guess=1;
c20guess=1/Param(ctrb).RMax;
end

for ctrb=1:K
  CoeffFileName=['Data/temp/War/c' ex(ctrb).casename '.mat'];
c10guess=1;
c20guess=1/Param(ctrb).RMax;
[sHist(:,ctrb),gHist(:,ctrb),u2btildHist(:,ctrb),RHist(:,ctrb),...
TauHist(:,ctrb),YHist(:,ctrb),TransHist(:,ctrb),btildHist(:,ctrb),...
c1Hist(:,ctrb),c2Hist(:,ctrb),l1Hist(:,ctrb),l2Hist(:,ctrb),...
IntHist(:,ctrb),IncomeFromAssets_Agent1Hist(:,ctrb),...
AfterTaxWageIncome_Agent1Hist(:,ctrb),AfterTaxWageIncome_Agent2Hist(:,ctrb),...
GShockDiffHist(:,ctrb),TransDiffHist(:,ctrb),LaborTaxAgent1DiffHist(:,ctrb),...
LaborTaxAgent2DiffHist(:,ctrb),DebtDiffHist(:,ctrb),GiniCoeffHist(:,ctrb)]...
=RunSimulations(CoeffFileName,0,c10guess,c20guess,NumSim,Param(ctrb),sHist0);
end

save( [Para.datapath 'SimDataParallelWar.mat'],'sHist',...
       'gHist','u2btildHist','RHist','TauHist','YHist','TransHist',...
       'btildHist','c1Hist','c2Hist','l1Hist','l2Hist','Para','IntHist',...
       'AfterTaxWageIncome_Agent1Hist','AfterTaxWageIncome_Agent2Hist',...
       'IncomeFromAssets_Agent1Hist','GShockDiffHist','TransDiffHist',...
       'LaborTaxAgent1DiffHist','LaborTaxAgent2DiffHist','DebtDiffHist',...
       'GiniCoeffHist')
    
   
   
   
%  % -- PLOT DIAGNOSTICS -----------------------------------------
close all
clear all
clc
SimTitle={'War'};
SimDataPath= 'Data/Calibration/War/SimDataParallelWar.mat';
SimPlotPath='Graphs/Calibration/War/';
mkdir(SimPlotPath)
SimTexPath='Tex/Calibration/War/';
mkdir(SimTexPath)
%PlotParallelSimulationsCommonShocks(SimDataPath,SimTexPath,SimPlotPath,SimTitle)


 Para.datapath=['Data/'];
 Para.StoreFileName=['cWar2.mat'];
 GetPlotsForFinalSolution(Para)
