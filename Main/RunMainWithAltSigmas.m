 clc
 clear all
 close all
 SetPath

 %{ 
This file solves the G-S economy with BGP preferences of the form
 psi.c^(1-sigma)/1-sigma+(1-psi)log[1-l] with following calibrations

1.] The ratio of productivities is 3.3 and the low productivity is
 normalized to 1
2.] psi is choosen to get a average FE of about 0.5
3.] pareto wts are such that the no-shock problem gives a tax rate of about
    20 percent
4.] Government expenditures are about 11 and 13 percent of the output
5.] beta =0.9
%}

% - XXXXXXX ANMOL - CURRENTLY IT USES A LEGACY METHOD FOR GETTING THE
% BASELINE PARAMETERS. NOW THAT WE HAVE A STEADY STATE CODE WE CAN USE THIS
% TO TARGET SOME AGREED MOMENTS IN OBSERVABLES 
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
opts = optimset('TolFun', 1e-15);
z=fsolve(@(z) GetCalibrationFrischElasticity (z,AvfFETarget,theta_1,theta_2,tau,g_Y,n1,n2), [1 1 ], opts);
gamma=z(1);
Y=z(2);


% BASELINE GOVERNMENT EXPENDITURE LEVELS
g=g_Y*Y;

% BASELINE PSI
psi=1/(1+gamma);
% BASELINE DISCOUNT FACTOR

beta=.9;

% BASELINE PARETO WTS
alpha_1=0.69;
alpha_2=1-alpha_1;
Para.n1=n1;
Para.n2=n2;
alpha_1=alpha_1*Para.n1;
alpha_2=alpha_2*Para.n2;

% BASELINE PROBABILITY MATRIX
NewPh=.5;
Para.P=[1-NewPh NewPh;1-NewPh NewPh];

% POPULATE THE PARA STRUC WITH THE BASELINE VALUES
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
Para.datapath=[rootDir sl 'Data/temp/'];
mkdir(Para.datapath)
casename='sigma';
Para.StoreFileName=['c' casename '.mat'];
CoeffFileName=[Para.datapath Para.StoreFileName];
 
 %  --- SOLVE THE BELLMAN EQUATION --------------------------------------
Para.Niter=200; % MAXIMUM NUMBER OF ITERATION


% flagSetRGrid,flagSetxGrid  TAKES TWO VALUES : 0 IF DEFAULT GRID OR 1 FOR USERDEFINED
% GRID

Para.flagSetRGrid=1; 
Para.flagSetxGrid=1;
Para.xMin=-2.5;
Para.xMax=2.5;

% EXPERIMENT 1 : SIGMA=1
casename='sigmaLow';
Para.StoreFileName=['c' casename '.mat'];
CoeffFileName=[Para.datapath Para.StoreFileName]; 
Para.sigma = 1;
Para.RMin=2.2;
Para.RMax=3.5;
MainBellman(Para) 


% EXPERIMENT 2 : SIGMA=2
casename='sigmaMed';
Para.StoreFileName=['c' casename '.mat'];
CoeffFileName=[Para.datapath Para.StoreFileName]; 
Para.sigma = 2;
Para.RMin=3.5;
Para.RMax=4.5;
MainBellman(Para) 

% EXPERIMENT 3 : SIGMA=3
casename='sigmaHigh';
Para.StoreFileName=['c' casename '.mat'];
CoeffFileName=[Para.datapath Para.StoreFileName];  
Para.sigma = 3;
Para.RMin=4.5;
Para.RMax=5.5;
MainBellman(Para) 



% SIMULATE THE ECONOMY
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
NumSim=60000;
rHist0 = rand(NumSim,1);
K=3;
ex(1).casename='sigmaLow'; 
ex(2).casename='sigmaMed';
ex(3).casename='sigmaHigh';




for ctrb=1:K
CoeffFileName=[rootDir sl 'Data/temp/c' ex(ctrb).casename '.mat'];
Sol=load(CoeffFileName);
Param(ctrb)=Sol.Para;
end

parfor ctrb=1:K
  CoeffFileName=['Data/temp/c' ex(ctrb).casename '.mat'];
c10guess=1;
c20guess=.5;

  
  [sHist(:,ctrb),gHist(:,ctrb),xHist(:,ctrb),RHist(:,ctrb),...
TauHist(:,ctrb),YHist(:,ctrb),TransHist(:,ctrb),btildHist(:,ctrb),...
c1Hist(:,ctrb),c2Hist(:,ctrb),l1Hist(:,ctrb),l2Hist(:,ctrb),...
IntHist(:,ctrb),IncomeFromAssets_Agent1Hist(:,ctrb),...
AfterTaxWageIncome_Agent1Hist(:,ctrb),AfterTaxWageIncome_Agent2Hist(:,ctrb),...
GShockDiffHist(:,ctrb),TransDiffHist(:,ctrb),LaborTaxAgent1DiffHist(:,ctrb),...
LaborTaxAgent2DiffHist(:,ctrb),DebtDiffHist(:,ctrb),GiniCoeffHist(:,ctrb)]...
=RunSimulations(CoeffFileName,0,c10guess,c20guess,NumSim,Param(ctrb),rHist0);
end

save([ rootDir sl 'Data/temp/SimDataParallelCommonShocks.mat'],'sHist',...
       'gHist','xHist','RHist','TauHist','YHist','TransHist',...
       'btildHist','c1Hist','c2Hist','l1Hist','l2Hist','Param','IntHist',...
       'AfterTaxWageIncome_Agent1Hist','AfterTaxWageIncome_Agent2Hist',...
       'IncomeFromAssets_Agent1Hist','GShockDiffHist','TransDiffHist',...
       'LaborTaxAgent1DiffHist','LaborTaxAgent2DiffHist','DebtDiffHist',...
       'GiniCoeffHist')
      
