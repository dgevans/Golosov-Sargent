function  [sHist,gHist,xHist,RHist,TauHist,YHist,TransHist,...
          btildHist,c1Hist,c2Hist,l1Hist,l2Hist,IntHist,...
          IncomeFromAssets_Agent1Hist,AfterTaxWageIncome_Agent1Hist,...
          AfterTaxWageIncome_Agent2Hist,GShockDiffHist,TransDiffHist,...
          LaborTaxAgent1DiffHist,LaborTaxAgent2DiffHist,DebtDiffHist,...
          GiniCoeffHist,xMeanHist,RmeanHist]=RunSimulationsAlt(CoeffFileName,btild0,c10guess,c20guess,NumSim,Para,rHist0)
% This function plots the similation for NumSim periods starting brom
% btild0 and using coeff from endIter. If existing draw of s-shocks are to
% be used..use the argument sHist0
if nargin==7
    flagUseExistingShocks='yes';
    disp('Using existing shocks')
end
close all;
olddatapath=Para.datapath;
oldtexpath=Para.texpath;
oldplotpath=Para.plotpath;
plotpath=oldplotpath;
datapath=olddatapath;

 load(CoeffFileName)
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
sigma=Para.sigma;
% SOLVE THE T-0 PROBLEM given btild(-1)
btild_1=btild0;
disp('Computed V...Now solving V0(btild_1) where btild_1 is')
disp(btild_1)
% c1 and c2 solve
options=optimset('Display','off');
[x,~,exitflagv0,~,~] = fminunc(@(x)  getValue0(x, btild_1,1,Para,c,V),[ c10guess c20guess],options);
if ~(exitflagv0==1)
[x,~,exitflagv0,~,~] = fminunc(@(x)  getValue0(x, btild_1,1,Para,c,V),[ 1 1/Para.RMax],options);
end

if ~(exitflagv0==1)
    disp('Optimization failed for V0 once ..trying with fmincon')
    opts = optimset('Algorithm', 'interior-point', 'Display','off', ...
        'GradObj','off','GradConstr','off',...
        'MaxIter',1000, ...
        'TolX', Para.ctol/10, 'TolFun', Para.ctol, 'TolCon', Para.ctol,'MaxTime',200);
    lb=[0.001 0.001];
    ub=[10 10];
    %[x,fval,exitflagv0,output,lambda]  =fmincon(@(x) getValue0(x, btild_1,1,Para,c,V),[ x ],[],[],[],[],lb,ub,[],opts);
    [x,~,exitflagv0,output,lambda]  =ktrlink(@(x) getValue0(x, btild_1,1,Para,c,V),[ c10guess c20guess],[],[],[],[],lb,ub,[],opts);
    
end
c10 = x(1);
c20 = x(2);
R0=c10/c20;
TotalResources=(c10*n1+c20*n2+g(1));
FF=R0*theta_2/theta_1;
DenL2=n1*theta_1*FF+theta_2*n2;
l20=(TotalResources-n1*theta_1+n1*theta_1*FF)/(DenL2);
l10= 1-FF*(1-l20);
BracketTerm=l20/(1-l20)-(l10/(1-l10))*R0;
xprime0=(((1-psi)/(psi))*BracketTerm+btild_1/(beta*psi)+R0-1)*psi;
btildprime0=xprime0/(c20^-1*psi) ;
Rprime0=c20^(-1)/c10^(-1);


% RUN SIMULATION
gHist=zeros(NumSim,1);
xHist=zeros(NumSim,1);
xMeanHist = zeros(NumSim,1);
btildHist=zeros(NumSim,1);
RHist=zeros(NumSim,1);
RmeanHist=zeros(NumSim,1);
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
xMeanHist(1) = xprime0;
ul20=(1-psi)/(1-l20);
ul10=(1-psi)/(1-l10);
uc20=psi/c20;
uc10=psi/c10;
c1Hist(1)=c10;
c2Hist(1)=c20;
l1Hist(1)=l10;
l2Hist(1)=l20;
btildHist(1)=btildprime0;
TauHist(1)=1-(ul20/(theta_2*uc20));
TransHist(1)=c20-l20*ul20/uc20;
RHist(1)=Rprime0;
RmeanHist(1) = Rprime0;
YHist(1)=n1*c10+n2*c20+g(1);
sHist(1)=1;
gHist(1)=g(sHist(1));
AfterTaxWageIncome_Agent1Hist(1)=l10*ul10/uc10;
AfterTaxWageIncome_Agent2Hist(1)=l20*ul10/uc20;
tic
for i=1:NumSim-1
    if mod(i,500)==0
        disp('Running Simulation, t=')
        disp(i)
        toc
        tic
    end
    % ------STATE (t) - x,R,s_ ------------------------------------------
    x=xHist(i);
    xmean = xMeanHist(i);
    R=RHist(i);
    Rmean = RmeanHist(i);
    s_=sHist(i);
    
    % ----SOLVE THE BELLMAN EQUATION  ------------------------------------
    [PolicyRulesInit]=GetInitialApproxPolicy([x R s_] ,domain,PolicyRulesStore);
    [PolicyRules, ~,exitflag,~]=CheckGradNAG(x,R,s_,c,V,PolicyRulesInit,Para,0);
    [PolicyRulesInitMean]=GetInitialApproxPolicy([xmean Rmean s_] ,domain,PolicyRulesStore);
    [PolicyRulesMean, ~,exitflagMean,~]=CheckGradNAG(xmean,Rmean,s_,c,V,PolicyRulesInitMean,Para);
    
    %---------------------------------------------------------------------------
    % GET THE POLICY RULES -
    %PolicyRules=[c1_1 c1_2 c2_1 c2_2 l1(1) l1(2) l2(1) l2(2) btildprime c2_1^(-1)/c1_1^(-1) c2_2^(-1)/c1_2^(-1) xprime(1) xprime(2)]
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
    xprime=PolicyRules(end-1:end);
    % btildprime - x'/u_c2
    btildprime=PolicyRules(9:10);
    
    
    xprimeMean=PolicyRulesMean(end-1:end);
    RprimeMean=PolicyRulesMean(end-3:end-2);
    
    % Int-Rates
    IntNum=psi/c2Hist(i); %marginal utility of consumption (Agent 2) in s_
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
   
    ExitFlag(i)=exitflag;
    
   
   
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
    RmeanHist(i+1) = Para.P(s_,:)*RprimeMean';
    xHist(i+1)=xprime(sHist(i+1)) ;
    xMeanHist(i+1) = Para.P(s_,:)*xprimeMean';
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
    %      xHist(n)=xprime(sHist(i+1)) ;
    %  n=n+1;
    %  end
    
end


end


% c1Hist(3:5)-btildHist(3:5)-AfterTaxWageIncome_Agent1Hist(3:5)-IncomeFromAssets_Agent1Hist(2:4)
% c2Hist(3:5)-AfterTaxWageIncome_Agent2Hist(3:5)
% YHist(3:5)-n1*c1Hist(3:5)-n2*c2Hist(3:5)-gHist(3:5)
% YHist(3:5)-n1*l1Hist(3:5)*theta_1-n2*l2Hist(3:5)*theta_2
% gHist(3:5)+n1*btildHist(3:5)+TransHist(3:5)-TauHist(3:5).*(n1*theta_1*l1Hist(3:5)+n2*theta_2*l2Hist(3:5))-n1*btildHist(2:4).*IntHist(2:4)