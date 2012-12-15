function  [sHist,gHist,u2btildHist,RHist, u2btildNLHist, RNLHist,TauHist,YHist,TransHist,...
          btildHist,c1Hist,c2Hist,l1Hist,l2Hist,IntHist,...
          IncomeFromAssets_Agent1Hist,AfterTaxWageIncome_Agent1Hist,...
          AfterTaxWageIncome_Agent2Hist,GShockDiffHist,TransDiffHist,...
          LaborTaxAgent1DiffHist,LaborTaxAgent2DiffHist,DebtDiffHist,...
          GiniCoeffHist]=RunSimulationsFromT1Linear(CoeffFileName,x0,R0,NumSim,Para,rHist0)
% This function plots the similation for NumSim periods starting brom
% btild0 and using coeff from endIter. If existing draw of s-shocks are to
% be used..use the argument sHist0
if nargin==6
    flagUseExistingShocks='yes';
    disp('Using existing shocks')
else
    flagUseExistingShocks='no';
    disp('yes')
end
close all;
olddatapath=Para.datapath;
oldtexpath=Para.texpath;
oldplotpath=Para.plotpath;
plotpath=oldplotpath;
datapath=olddatapath;

load(CoeffFileName);

[A XSS] = LinearApproximation(Para);
xSS = XSS(11);
RSS = XSS(9);

disp('Govt Exp')
g=Para.g
n1=Para.n1;
n2=Para.n2;
alpha_1=Para.alpha_1;
alpha_2=Para.alpha_2;
%disp('Govt Exp')
%g=Para.g
theta_1=Para.theta_1;
theta_2=Para.theta_2;
psi=Para.psi;
beta=Para.beta;
sigma = Para.sigma;

% RUN SIMULATION
gHist=zeros(NumSim,1);
u2btildHist=zeros(NumSim,1);
u2btildNLHist=zeros(NumSim,1);
btildHist=zeros(NumSim,1);
RHist=zeros(NumSim,1);
RNLHist=zeros(NumSim,1);
TauHist=zeros(NumSim,1);
YHist=zeros(NumSim,1);
TransHist=zeros(NumSim,1);
GMul=zeros(NumSim,1);
c1Hist=zeros(NumSim,1);
c2Hist=zeros(NumSim,1);
l1Hist=zeros(NumSim,1);
l2Hist=zeros(NumSim,1);
sHist=zeros(NumSim,1);
GiniCoeffHist=zeros(NumSim,1);
IntHist=zeros(NumSim-1,1);
IncomeFromAssets_Agent1Hist=zeros(NumSim-1,1);
AfterTaxWageIncome_Agent1Hist=zeros(NumSim,1);
AfterTaxWageIncome_Agent2Hist=zeros(NumSim,1);
GShockDiffHist=zeros(NumSim-1,1);
TransDiffHist=zeros(NumSim-1,1);
LaborTaxAgent1DiffHist=zeros(NumSim-1,1);
LaborTaxAgent2DiffHist=zeros(NumSim-1,1);
DebtDiffHist=zeros(NumSim-1,1);
sHist(1)=1;


% INITIALIZE - t=0
u2btildHist(1)=x0;
RHist(1)=R0;
u2btildNLHist(1)=x0;
RNLHist(1)=R0;
tic
for i=1:NumSim-1
    if mod(i,500)==0
        disp('Running Simulation, t=')
        disp(i)
        toc
        tic
    end
    % ------STATE (t) - x,R,s_ ------------------------------------------
    u2btild=u2btildHist(i);
    R=RHist(i);
    u2btildNL=u2btildNLHist(i);
    RNL=RNLHist(i);
    s_=sHist(i);
    
    % ----Use linear Approximation to get policies  ------------------------------------
    yhat = [R; u2btild] - [RSS; xSS];
    PolicyRules = (A*yhat)'+XSS;
    [PolicyRulesInit]=GetInitialApproxPolicy([u2btildNL RNL s_] ,x_state,PolicyRulesStore);
    [PolicyRulesNL, ~,exitflag,~]=CheckGradNAG(u2btildNL,RNL,s_,c,V,PolicyRulesInit,Para,0);
    
    %---------------------------------------------------------------------------
    % GET THE POLICY RULES -
    %PolicyRules=[c1_1 c1_2 c2_1 c2_2 l1(1) l1(2) l2(1) l2(2) btildprime c2_1^(-1)/c1_1^(-1) c2_2^(-1)/c1_2^(-1) u2btildprime(1) u2btildprime(2)]
    c1=PolicyRules(1:2);
    c2=PolicyRules(3:4);
    l1=PolicyRules(5:6);
    l2=PolicyRules(7:8);
    ul2=(1-psi)./(1-l2);
    uc2=psi.*c2.^(-sigma);
    ul1=(1-psi)./(1-l1);
    uc1=psi.*c1.^(-sigma);
    Rprime=PolicyRules(9:10);
    % x' - u_c_2* btildprime
    u2btildprime=PolicyRules(11:12);
    RprimeNL=PolicyRulesNL(end-3:end-2);
    % x' - u_c_2* btildprime
    u2btildprimeNL=PolicyRulesNL(end-1:end);
    % btildprime - x'/u_c2
    btildprime=u2btildprime./uc2;
    
    % Int-Rates
    IntNum=psi.*c2Hist(i).^(-sigma); %marginal utility of consumption (Agent 2) in s_
    DenNum=beta*sum(Para.P(sHist(i),:).*uc2); % expected marginal utility of consumption (Agent 2)
    IntHist(i)=IntNum/DenNum;
    
   
    % TAU - From the WAGE optimality of Agent 2
    Tau=1-(ul2./(theta_2.*uc2));
    
    % OUTPUT
    y(1)=c1(1)*n1+c2(1)*n2+g(1);
    y(2)=c1(2)*n1+c2(2)*n2+g(2);
    
    % TRANSFERS
    % These are transfers computed on the assumption that Agent 2 cannot
    % borrow and lend. The transfers are the difference between his
    % consumption and after tax earning (l . U_l/U_c)
    Trans=c2-l2.*ul2./uc2;
    
    % G MULTIPLIER - Computed using (yh-yl)/(gh-gl)
    GMul(i)=(y(2)-y(1))/(g(2)-g(1));
    
     % Income
    AfterTaxWageIncome_Agent2=l2.*ul2./uc2+Trans;
    AfterTaxWageIncome_Agent1=l1.*ul1./uc1+Trans;
    % Gini Coeff
    GiniCoeff=(AfterTaxWageIncome_Agent2 +2*AfterTaxWageIncome_Agent1)./(AfterTaxWageIncome_Agent2+AfterTaxWageIncome_Agent1)-3/2;
   
    
   
   
    % DRAW THE s' ~ P(s,:) if flagUseExistingShocks is set to no
    if strcmpi(flagUseExistingShocks,'yes')
    if rHist0(i+1) < Para.P(sHist(i),1)
        sHist(i+1)=1;
    else
        
        sHist(i+1)=2;
    end
    else
    if rand < Para.P(sHist(i),1)
        sHist(i+1)=1;
    else
        
        sHist(i+1)=2;
    end
    end
    % UPDATE THE SIMULATION HISTORY
    
    RHist(i+1)=Rprime(sHist(i+1));
    u2btildHist(i+1)=u2btildprime(sHist(i+1)) ;
    RNLHist(i+1)=RprimeNL(sHist(i+1));
    u2btildNLHist(i+1)=u2btildprimeNL(sHist(i+1)) ;
    btildHist(i+1)=btildprime(sHist(i+1)) ;
    TauHist(i+1)=Tau(sHist(i+1));
    YHist(i+1)=y(sHist(i+1));
    TransHist(i+1)=Trans(sHist(i+1));
    c1Hist(i+1)=c1(sHist(i+1));
    c2Hist(i+1)=c2(sHist(i+1));
    l1Hist(i+1)=l1(sHist(i+1));
    l2Hist(i+1)=l2(sHist(i+1));
    IncomeFromAssets_Agent1Hist(i)=-btildHist(i).*(IntHist(i)-1);
    AfterTaxWageIncome_Agent1Hist(i+1)=AfterTaxWageIncome_Agent1(sHist(i+1));
    AfterTaxWageIncome_Agent2Hist(i+1)=AfterTaxWageIncome_Agent2(sHist(i+1));
    gHist(i+1)=g(sHist(i+1));
     % Diff in GBC
    % diff in g_shock
    GShockDiffHist(i)=g(2)-g(1);
    % diff in trasnfers
    TransDiffHist(i)=(Trans(2)-Trans(1));
    % diff in labortax agent1
    LaborTaxAgent1DiffHist(i)=theta_1*l1(2)*Tau(2)*n1 - theta_1*l1(1)*Tau(1)*n1;
    % diff in labortax agent2
    LaborTaxAgent2DiffHist(i)=theta_2*l2(2)*Tau(2)*n2 - theta_2*l2(1)*Tau(1)*n2;
      % diff in borrowing
    DebtDiffHist(i)=n1*(btildprime(2)-btildprime(1));
    GiniCoeffHist(i+1)=GiniCoeff(sHist(i+1));
    %  if exitflag==1
    %      RHist(n)=Rprime(sHist(i+1));
    %      u2btildHist(n)=u2btildprime(sHist(i+1)) ;
    %  n=n+1;
    %  end
    
end


end


% c1Hist(3:5)-btildHist(3:5)-AfterTaxWageIncome_Agent1Hist(3:5)-IncomeFromAssets_Agent1Hist(2:4)
% c2Hist(3:5)-AfterTaxWageIncome_Agent2Hist(3:5)
% YHist(3:5)-n1*c1Hist(3:5)-n2*c2Hist(3:5)-gHist(3:5)
% YHist(3:5)-n1*l1Hist(3:5)*theta_1-n2*l2Hist(3:5)*theta_2
% gHist(3:5)+n1*btildHist(3:5)+TransHist(3:5)-TauHist(3:5).*(n1*theta_1*l1Hist(3:5)+n2*theta_2*l2Hist(3:5))-n1*btildHist(2:4).*IntHist(2:4)