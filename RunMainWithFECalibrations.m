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
alpha_1=0.69;
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
Para.theta_1=theta_1;
Para.theta_2=theta_2;
Para.btild_1=0;
Para.alpha_1=alpha_1;
Para.alpha_2=alpha_2;
Para.datapath=['Data/temp/'];
mkdir(Para.datapath)
casename='FE_High';
Para.StoreFileName=['c' casename '.mat'];
CoeffFileName=[Para.datapath Para.StoreFileName];
 
 %  --- SOLVE THE BELLMAN EQUATION --------------------------------------
 % test run 
 Para.Niter=250;
 Para.sigma = 1;
RGrid.RMin=2.2;
RGrid.RMax=3.2;
NewPh=.5;
Para.P=[1-NewPh NewPh;1-NewPh NewPh];
MainBellmanAltInit(Para,RGrid) 


% --- Med alpha ---------------------------------------------------------

alpha_1=0.5;
alpha_2=1-alpha_1;
alpha_1=alpha_1*Para.n1;
alpha_2=alpha_2*Para.n2;
Para.alpha_1=alpha_1;
Para.alpha_2=alpha_2;
casename='FE_Med';
Para.StoreFileName=['c' casename '.mat'];
RGrid.RMin=2.2;
RGrid.RMax=2.8;
%MainBellman(Para,RGrid)



% --- Low alpha ---------------------------------------------------------

alpha_1=0.25;
alpha_2=1-alpha_1;
alpha_1=alpha_1*Para.n1;
alpha_2=alpha_2*Para.n2;
Para.alpha_1=alpha_1;
Para.alpha_2=alpha_2;
casename='FE_Low';
Para.StoreFileName=['c' casename '.mat'];
RGrid.RMin=2.2;
RGrid.RMax=2.5;
%MainBellman(Para,RGrid)



 %% Set the Parallel Config
err=[];
try
    matlabpool('size')
catch err
end
if isempty(err)
    
    
    if(matlabpool('size') > 0)
        matlabpool close
    end
    
    matlabpool open local;
    
end

%-- Simulate the MODEL -------------------------------------------------
NumSim=10000;
rHist0 = rand(NumSim,1);


K=3;

ex(1).casename='PhMed'; % benchmark calibrations high alpha1
ex(2).casename='PhHigh'; % benchmark calibrations with medium alpha1
ex(3).casename='PhHighHigh'; % benchmark calibrations high alpha1



for ctrb=1:K
CoeffFileName=['Data/temp/c' ex(ctrb).casename '.mat'];
Sol=load(CoeffFileName);
Param(ctrb)=Sol.Para;
end

parfor ctrb=1:K
  CoeffFileName=['Data/temp/c' ex(ctrb).casename '.mat'];
c10guess=1;
c20guess=1/Param(ctrb).RMax;

  
  [sHist(:,ctrb),gHist(:,ctrb),u2btildHist(:,ctrb),RHist(:,ctrb),...
TauHist(:,ctrb),YHist(:,ctrb),TransHist(:,ctrb),btildHist(:,ctrb),...
c1Hist(:,ctrb),c2Hist(:,ctrb),l1Hist(:,ctrb),l2Hist(:,ctrb),...
IntHist(:,ctrb),IncomeFromAssets_Agent1Hist(:,ctrb),...
AfterTaxWageIncome_Agent1Hist(:,ctrb),AfterTaxWageIncome_Agent2Hist(:,ctrb),...
GShockDiffHist(:,ctrb),TransDiffHist(:,ctrb),LaborTaxAgent1DiffHist(:,ctrb),...
LaborTaxAgent2DiffHist(:,ctrb),DebtDiffHist(:,ctrb),GiniCoeffHist(:,ctrb),...
u2btildMeanHist(:,ctrb),RMeanHist(:,ctrb)]...
=RunSimulationsAlt(CoeffFileName,0,c10guess,c20guess,NumSim,Param(ctrb),rHist0);
end

save( [Para.datapath 'SimDataParallelPertPAlt.mat'],'sHist',...
       'gHist','u2btildHist','RHist','TauHist','YHist','TransHist',...
       'btildHist','c1Hist','c2Hist','l1Hist','l2Hist','Para','IntHist',...
       'AfterTaxWageIncome_Agent1Hist','AfterTaxWageIncome_Agent2Hist',...
       'IncomeFromAssets_Agent1Hist','GShockDiffHist','TransDiffHist',...
       'LaborTaxAgent1DiffHist','LaborTaxAgent2DiffHist','DebtDiffHist',...
       'GiniCoeffHist','u2btildMeanHist','RMeanHist')

   
   
   
%  % -- PLOT DIAGNOSTICS -----------------------------------------
close all
clear all
clc
SimTitle{1}='$\alpha_1=0.69$';
SimTitle{2}='$\alpha_1=0.50$';
SimTitle{3}='$\alpha_1=0.25$';
SimDataPath= 'Data/Calibration/SimDataParallelCommonShocks.mat';
SimPlotPath='Graphs/Calibration/';
mkdir(SimPlotPath)
SimTexPath='Tex/Calibration/';
mkdir(SimTexPath)
PlotParallelSimulationsCommonShocks(SimDataPath,SimTexPath,SimPlotPath,SimTitle)


 Para.datapath=['Data/temp/'];
 %Para.StoreFileName=['c' ex(2).casename '.mat'];
 Para.StoreFileName=['cPhHigh.mat'];
 
 GetPlotsForFinalSolution(Para)
 for i=1:3
SimDataPath= 'Data/Calibration/SimDataParallelCommonShocks.mat';
load(SimDataPath)
Domain.xBounds=[min(u2btildHist(:,i)) max(u2btildHist(:,i))];
Domain.RBounds=[min(RHist(:,i)) max(RHist(:,i))];
 Para.datapath=['Data/Calibration/'];
 Para.StoreFileName=['c' ex(i).casename '.mat'];
 plotpath=['Graphs/Calibration/' ex(i).casename '/']    
GetPlotsForFinalSolution(Para,plotpath)
 end