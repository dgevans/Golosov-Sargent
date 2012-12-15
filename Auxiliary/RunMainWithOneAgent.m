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
g_l_y=.12; % g low
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
 Para.Niter=150;
RGrid.RMin=3;
RGrid.RMax=5;
Para.flagSetu2BtildGrid=1;
Para.u2btildMin=-3;
Para.u2btildMax=3;
matlabpool close force local
Indx=MainBellman(Para,RGrid) 
Indx=150
matlabpool close force local
Para.StoreFileName=['c_' num2str(Indx) '.mat'];
Para.flagPlot2PeriodDrifts=0
GetPlotsForFinalSolution(Para)
load([Para.datapath Para.StoreFileName])
xState=fsolve(@(x) GetCrossingPoints(x,1,c,V,PolicyRulesStore,x_state,Para),[mean(Para.u2bdiffGrid) mean(Para.RGrid)])
u2btild=xState(1);
 R=xState(2);
 [xInit]=GetInitialApproxPolicy([u2btild R 1],x_state,PolicyRulesStore);
            %xInit=PolicyRulesStore(ctr,:);
            [PolicyRules, V_new,exitflag,~]=CheckGradNAG(u2btild,R,1,c,V,xInit',Para,0);
        psi=Para.psi;
        beta=Para.beta;
        g=Para.g;
        theta_1=Para.theta_1;
        n1=Para.n1;
        n2=Para.n2;
        c1=PolicyRules(1:2);
    c2=PolicyRules(3:4);
    l1=PolicyRules(5:6);
    l2=PolicyRules(7:8);
    ul2=(1-psi)./(1-l2);
    uc2=psi./c2;
    ul1=(1-psi)./(1-l1);
    uc1=psi./c1;
    Rprime=PolicyRules(end-3:end-2);
    % x' - u_c_2* btildprime
    u2btildprime=PolicyRules(end-1:end);
    % btildprime - x'/u_c2
    btildprime=PolicyRules(9:10);
    
    % Int-Rates
    
    % TAU - From the WAGE optimality of Agent 2
    Tau=1-(ul1./(theta_1.*uc1));
    
    % OUTPUT
    y(1)=c1(1)*n1+c2(1)*n2+g(1);
    y(2)=c1(2)*n1+c2(2)*n2+g(2);
    
    % TRANSFERS
    
    u2btild-btildprime
    (theta_1*l1)./(u2btild./(beta*sum(Para.P(1,:).*c2.^-1))-btildprime)
    % These are transfers computed on the assumption that Agent 2 cannot
    % borrow and lend. The transfers are the difference between his
    % consumption and after tax earning (l . U_l/U_c)
    Trans=c2;
    
    
     % Income
    AfterTaxWageIncome_Agent2=Trans;
    AfterTaxWageIncome_Agent1=l1.*ul1./uc1+Trans;
    % Gini Coeff
    GiniCoeff=(AfterTaxWageIncome_Agent2 +2*AfterTaxWageIncome_Agent1)./(AfterTaxWageIncome_Agent2+AfterTaxWageIncome_Agent1)-3/2;
 
%-- Simulate the MODEL -------------------------------------------------
NumSim=25000;
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

   
 