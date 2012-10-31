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
theta_1=2; % theta high
theta_2=0.001;  % theta low
g_l_y=.125; % g low
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
alpha_1=0.7;
alpha_2=1-alpha_1;
theta_1=2; % theta high
theta_2=0;  % theta low
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
casename='OneAgent';
Para.StoreFileName=['c' casename '.mat'];
CoeffFileName=[Para.datapath Para.StoreFileName];
 
 %  --- SOLVE THE BELLMAN EQUATION --------------------------------------
 % test run 
 Para.Niter=250;
RGrid.RMin=1;
RGrid.RMax=1.5;
Para.flagSetu2BtildGrid=1;
Para.u2btildMin=1;
Para.u2btildMax=2;
matlabpool close force local
Indx=MainBellman(Para,RGrid) 
Indx=5
matlabpool close force local
Para.StoreFileName=['c_' num2str(Indx) '.mat'];
Para.flagPlot2PeriodDrifts=0
GetPlotsForFinalSolution(Para)

%-- Simulate the MODEL -------------------------------------------------
NumSim=100;
rHist0 = rand(NumSim,1);


K=1;

ex(1).casename='OneAgent'; % benchmark calibrations high alpha1



for ctrb=1:K
CoeffFileName=['Data/temp/c' ex(ctrb).casename '.mat'];
Sol=load(CoeffFileName);
Param(ctrb)=Sol.Para;
end

for ctrb=1:K
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
]...
=RunSimulations(CoeffFileName,(Para.u2btildMin+Para.u2btildMax)/2,c10guess,c20guess,NumSim,Param(ctrb),rHist0);
end

save( [Para.datapath 'SimDataParallelPertPAlt.mat'],'sHist',...
       'gHist','u2btildHist','RHist','TauHist','YHist','TransHist',...
       'btildHist','c1Hist','c2Hist','l1Hist','l2Hist','Para','IntHist',...
       'AfterTaxWageIncome_Agent1Hist','AfterTaxWageIncome_Agent2Hist',...
       'IncomeFromAssets_Agent1Hist','GShockDiffHist','TransDiffHist',...
       'LaborTaxAgent1DiffHist','LaborTaxAgent2DiffHist','DebtDiffHist',...
       'GiniCoeffHist')

   
 