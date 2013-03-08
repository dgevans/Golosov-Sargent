function  [SimData]=RunSimulations(CoeffFileName,btild0,c10guess,c20guess,NumSim,Para,rHist0)
      
% THIS FUCNTION COMPUTES THE SIMLATION USING THE USER GIVEN SEED
if nargin==7
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
tempSavePath=Para.saveSimPath;
 load(CoeffFileName)
disp('Govt Exp')
g=Para.g;
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
sigma=Para.sigma;
P=Para.P;
S=length(P(1,:));
CumP=[];

for s_ = 1:S
    CumP(s_,1)=0;
for s=2:S
    CumP(s_,s)=CumP(s_,s-1)+P(s_,s);
end
CumP(s_,S+1)=1;
end


% SOLVE THE T-0 PROBLEM given btild(-1)
btild_1=btild0;
s_=1;

disp('Computed V...Now solving V0(btild_1) where btild_1 is')
disp(btild_1)
% c1 and c2 solve
options=optimset('Display','off');
[x,~,exitflagv0,~,~] = fminunc(@(x)  getValue0(x, btild_1,s_,Para,c,V),[ c10guess c20guess],options);
if ~(exitflagv0==1)
[x,~,exitflagv0,~,~] = fminunc(@(x)  getValue0(x, btild_1,s_,Para,c,V),[ 1 1/Para.RMax],options);
end

if ~(exitflagv0==1)
    disp('Optimization failed for V0 once ..trying with fmincon')
    opts = optimset('Algorithm', 'interior-point', 'Display','off', ...
        'GradObj','off','GradConstr','off',...
        'MaxIter',1000, ...
        'TolX', Para.ctol/10, 'TolFun', Para.ctol, 'TolCon', Para.ctol,'MaxTime',200);
    lb=[0.001 0.001];
    ub=[10 10];
    [x,fval,exitflagv0,output,lambda]  =fmincon(@(x) getValue0(x, btild_1,s_,Para,c,V),[ x ],[],[],[],[],lb,ub,[],opts);
    %[x,~,exitflagv0,output,lambda]  =ktrlink(@(x) getValue0(x, btild_1,1,Para,c,V),[ c10guess c20guess],[],[],[],[],lb,ub,[],opts);
    
end
c10 = x(1);
c20 = x(2);
R0=(c10/c20)^(sigma);
TotalResources=(c10*n1+c20*n2+g(s_));
DenL2=theta_2*R0*n1+theta_2*n2;
l20=(TotalResources-theta_1*n1+ theta_2*n1*R0)/(DenL2);
l10= 1-(1-l20)*theta_2/theta_1*R0;
xprime0=-(c20-c10)*(psi*c20^(-sigma))-((l10/(1-l10))*R0-l20/(1-l20))*(1-psi)+btild_1*psi*c20^(-sigma);
Rprime0=c20^(-sigma)/c10^(-sigma);

% RUN SIMULATION
gHist=zeros(NumSim,1);
xHist=zeros(NumSim,1);
btildHist=zeros(NumSim,1);
RHist=zeros(NumSim,1);
btildHist(1)=btild_1;
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


% INITIALIZE - t=0
xHist(1)=xprime0;
ul20=(1-psi)/(1-l20);
ul10=(1-psi)/(1-l10);
uc20=psi/(c20^(sigma));
uc10=psi/(c10^(sigma));
c1Hist(1)=c10;
c2Hist(1)=c20;
l1Hist(1)=l10;
l2Hist(1)=l20;
btildHist(1)=xprime0/uc20;
TauHist(1)=1-(ul10./(theta_1*uc10));
TransHist(1)=c20-l20*ul20/uc20;
RHist(1)=Rprime0;
YHist(1)=n1*c10+n2*c20+g(s_);
sHist(1)=s_;
gHist(1)=g(sHist(1));
AfterTaxWageIncome_Agent1Hist(1)=l10*ul10/uc10;
AfterTaxWageIncome_Agent2Hist(1)=l20*ul10/uc20;
tic
for i=1:NumSim-1
    % ------STATE (t) - x,R,s_ ------------------------------------------
    x=xHist(i);
    R=RHist(i);
    s_=sHist(i);
    
    % ----SOLVE THE BELLMAN EQUATION  ------------------------------------
    [PolicyRulesInit]=GetInitialApproxPolicy([x R s_] ,domain,PolicyRulesStore);
    [PolicyRules, ~,exitflag,~]=CheckGradNAG(x,R,s_,c,V,PolicyRulesInit,Para);
    
    %---------------------------------------------------------------------------
    % GET THE POLICY RULES -
    %PolicyRules=[c1_1 c1_2 c2_1 c2_2 l1(1) l1(2) l2(1) l2(2) btildprime c2_1^(-1)/c1_1^(-1) c2_2^(-1)/c1_2^(-1) xprime(1) xprime(2)]
    c1=PolicyRules(1:S);
    c2=PolicyRules(S+1:2*S);
    l1=PolicyRules(2*S+1:3*S);
    l2=PolicyRules(3*S+1:4*S);
    ul2=(1-psi)./(1-l2);
    uc2=psi./(c2.^(sigma));
    ul1=(1-psi)./(1-l1);
    uc1=psi./(c1.^(sigma));
    Rprime=PolicyRules(end-2*S+1:end-S);
    % x' - u_c_2* btildprime
    xprime=PolicyRules(end-S+1:end);
    % btildprime - x'/u_c2
    btildprime=PolicyRules(end-3*S+1:end-2*S);
    
    % Int-Rates
    IntNum=psi/(c2Hist(i)^(sigma)); %marginal utility of consumption (Agent 2) in s_
    DenNum=beta*sum(Para.P(sHist(i),:).*uc2); % expected marginal utility of consumption (Agent 2)
    IntHist(i)=IntNum/DenNum;
    
   
    % TAU - From the WAGE optimality of Agent 2
    Tau=1-(ul1./(theta_1.*uc1));
    
    % OUTPUT
    
    y=c1*n1+c2*n2+g;
    
    
    % TRANSFERS
    % These are transfers computed on the assumption that Agent 2 cannot
    % borrow and lend. The transfers are the difference between his
    % consumption and after tax earning (l . U_l/U_c)
    Trans=c2-l2.*ul2./uc2;
    
    % G MULTIPLIER - Computed using (yh-yl)/(gh-gl)
    GMul(i)=(y(S)-y(1))/(g(S)-g(1));
    
     % Income
    AfterTaxWageIncome_Agent2=l2.*ul2./uc2+Trans;
    AfterTaxWageIncome_Agent1=l1.*ul1./uc1+Trans;
    % Gini Coeff
    GiniCoeff=(AfterTaxWageIncome_Agent2 +2*AfterTaxWageIncome_Agent1)./(AfterTaxWageIncome_Agent2+AfterTaxWageIncome_Agent1)-3/2;
   
    ExitFlag(i)=exitflag;
    
   
   
    % DRAW THE s' ~ P(s,:) if flagUseExistingShocks is set to no
    if strcmpi(flagUseExistingShocks,'yes')
        sHist(i+1)=sum(~(CumP(sHist(i),:)-rHist0(i+1)>0));

    else
            sHist(i+1)=sum(~(CumP(sHist(i),:)-rand>0));
    end
    % UPDATE THE SIMULATION HISTORY
    
    RHist(i+1)=Rprime(sHist(i+1));
    xHist(i+1)=xprime(sHist(i+1)) ;
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
    GShockDiffHist(i)=g(S)-g(1);
    % diff in trasnfers
    TransDiffHist(i)=(Trans(S)-Trans(1));
    % diff in labortax agent1
    LaborTaxAgent1DiffHist(i)=theta_1*l1(S)*Tau(S)*n1 - theta_1*l1(1)*Tau(1)*n1;
    % diff in labortax agent2
    LaborTaxAgent2DiffHist(i)=theta_2*l2(S)*Tau(S)*n2 - theta_2*l2(1)*Tau(1)*n2;
      % diff in borrowing
    DebtDiffHist(i)=n1*(btildprime(S)-btildprime(1));
    GiniCoeffHist(i+1)=GiniCoeff(sHist(i+1)); 
   
    
     if mod(i,1000)==0 || i==NumSim-1
        disp('Running Simulation, t=')
        disp(i)
        toc
        tic
        SimData.sHist=sHist;
SimData.gHist=gHist;
SimData.xHist=xHist;
SimData.RHist=RHist;
SimData.TauHist=TauHist;
SimData.YHist=YHist;
SimData.TransHist=TransHist;
SimData.btildHist=btildHist;
SimData.c1Hist=c1Hist;
SimData.c2Hist=c2Hist;
SimData.l1Hist=l1Hist;
SimData.l2Hist=l2Hist;
SimData.IntHist=IntHist;
SimData.IncomeFromAssets_Agent1Hist=IncomeFromAssets_Agent1Hist;
SimData.AfterTaxWageIncome_Agent1Hist=AfterTaxWageIncome_Agent1Hist;
SimData.AfterTaxWageIncome_Agent2Hist=AfterTaxWageIncome_Agent2Hist;
SimData.GShockDiffHist=GShockDiffHist;
SimData.TransDiffHist=TransDiffHist;
SimData.LaborTaxAgent1DiffHist=LaborTaxAgent1DiffHist;
SimData.LaborTaxAgent2DiffHist=LaborTaxAgent2DiffHist;
SimData.DebtDiffHist=DebtDiffHist;
SimData.GiniCoeffHist=GiniCoeffHist;
save(tempSavePath,'SimData')
    end
   
end

end

% BUDGET CONSTRAINTS
% c1Hist(3:5)-btildHist(3:5)-AfterTaxWageIncome_Agent1Hist(3:5)-IncomeFromAssets_Agent1Hist(2:4)
% c2Hist(3:5)-AfterTaxWageIncome_Agent2Hist(3:5)
% YHist(3:5)-n1*c1Hist(3:5)-n2*c2Hist(3:5)-gHist(3:5)
% YHist(3:5)-n1*l1Hist(3:5)*theta_1-n2*l2Hist(3:5)*theta_2
% gHist(3:5)+n1*btildHist(3:5)+TransHist(3:5)-TauHist(3:5).*(n1*theta_1*l1Hist(3:5)+n2*theta_2*l2Hist(3:5))-n1*btildHist(2:4).*IntHist(2:4)