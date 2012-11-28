 clc
 clear all
 close all


%-- Simulate the MODEL -------------------------------------------------
NumSim=25000;
rHist0=(rand(NumSim,1));

x0=3;
R0=6.2;
K=1;

ex(1).casename='HighAlpha'; % benchmark calibrations high alpha1



for ctrb=1:K
CoeffFileName=['Data/Calibration/c' ex(ctrb).casename '.mat'];
Sol=load(CoeffFileName);
Param(ctrb)=Sol.Para;
end

for ctrb=1:K
  CoeffFileName=['Data/Calibration/c' ex(ctrb).casename '.mat'];
c10guess=1;
c20guess=1/Param(ctrb).RMax;
[sHist(:,ctrb),gHist(:,ctrb),u2btildHist(:,ctrb),RHist(:,ctrb),...
TauHist(:,ctrb),YHist(:,ctrb),TransHist(:,ctrb),btildHist(:,ctrb),...
c1Hist(:,ctrb),c2Hist(:,ctrb),l1Hist(:,ctrb),l2Hist(:,ctrb),...
IntHist(:,ctrb),IncomeFromAssets_Agent1Hist(:,ctrb),...
AfterTaxWageIncome_Agent1Hist(:,ctrb),AfterTaxWageIncome_Agent2Hist(:,ctrb),...
GShockDiffHist(:,ctrb),TransDiffHist(:,ctrb),LaborTaxAgent1DiffHist(:,ctrb),...
LaborTaxAgent2DiffHist(:,ctrb),DebtDiffHist(:,ctrb),GiniCoeffHist(:,ctrb)]...
=RunSimulationsFromT1(CoeffFileName,x0,R0,NumSim,Param(ctrb),rHist0);
end

save( ['Data/Calibration/SimulationHighAlpha.mat'],'sHist',...
       'gHist','u2btildHist','RHist','TauHist','YHist','TransHist',...
       'btildHist','c1Hist','c2Hist','l1Hist','l2Hist','IntHist',...
       'AfterTaxWageIncome_Agent1Hist','AfterTaxWageIncome_Agent2Hist',...
       'IncomeFromAssets_Agent1Hist','GShockDiffHist','TransDiffHist',...
       'LaborTaxAgent1DiffHist','LaborTaxAgent2DiffHist','DebtDiffHist',...
       'GiniCoeffHist')

   
   
 