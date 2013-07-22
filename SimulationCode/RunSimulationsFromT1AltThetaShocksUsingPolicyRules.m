function  [SimData]=RunSimulationsFromT1AltThetaShocksUsingPolicyRules(PolicyFunctions,InitialConditions,NumSim,Para,rHist0)
% This function plots the similation for NumSim periods starting brom
% btild0 and using coeff from endIter. If existing draw of s-shocks are to
% be used..use the argument sHist0


% grab the policy functions
c1=PolicyFunctions.c1;
c2=PolicyFunctions.c2;
l1=PolicyFunctions.l1;
l2=PolicyFunctions.l2;
xPrime=PolicyFunctions.xPrime;
RPrime=PolicyFunctions.RPrime;



close all;

 disp('theta Exp')

 % grab the parameters 
 
 g=Para.g;
n1=Para.n1;
n2=Para.n2;
%disp('Govt Exp')
%g=Para.g
theta_1=Para.theta_1;
theta_2=Para.theta_2;
psi=Para.psi;
beta=Para.beta;
sigma=Para.sigma;
P=Para.P;
S=length(P(1,:));
CumP=zeros(S,S);

for s_ = 1:S
    for s=1:S
    CumP(s_,s)=CumP(s_,s)+P(s_,s);
    end
end


% RUN SIMULATION
Theta_1Hist=zeros(NumSim,1);
Theta_2Hist=zeros(NumSim,1);
xHist=zeros(NumSim,1);
btildHist=zeros(NumSim,1);
RHist=zeros(NumSim,1);
TauHist=zeros(NumSim,1);
YHist=zeros(NumSim,1);
TransHist=zeros(NumSim,1);
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




% INITIALIZE - t=0
sHist(1)=InitialConditions.s0;
xHist(1)=InitialConditions.x0;
RHist(1)=InitialConditions.R0;


tic

for i=1:NumSim-1
    % ------STATE (t) - x,R,s_ ------------------------------------------
    x=xHist(i);
    R=RHist(i);
    s_=sHist(i);
    
    for s=1:S
    if rHist0(i)<CumP(s_,s)
        break;
    end
    end
    sHist(i+1)=s;
    %sHist(i+1)=discretesample(P(s_,:),1);
           
    for s=1:S
    C1(s)=c1{s_,s}(x,R);
    C2(s)=c2{s_,s}(x,R);
    L1(s)=l1{s_,s}(x,R);
    L2(s)=l2{s_,s}(x,R);
  RRPrime(s)=max(min(RPrime{s_,s}(x,R),Para.RMax),Para.RMin);
    %RRPrime(s)=RPrime{s_,s}(x,R);
    xxPrime(s)=xPrime{s_,s}(x,R);
    end
    
    
    ul2=(1-psi)./(1-L2);
    uc2=psi./(C2.^(sigma));
    ul1=(1-psi)./(1-L1);
    uc1=psi./(C1.^(sigma));
    btildprime=xxPrime./uc2;
    
    % Int-Rates
    IntNum=psi/(c2Hist(i)^(sigma)); %marginal utility of consumption (Agent 2) in s_
    DenNum=beta(s_)*sum(Para.P(s_,:).*uc2); % expected marginal utility of consumption (Agent 2)
    IntHist(i)=IntNum./DenNum;
    
   
    % TAU - From the WAGE optimality of Agent 1
    Tau=1-(ul1./(theta_1.*uc1));
    
    % OUTPUT
    
    y=C1*n1+C2*n2+g;
    
    
    % TRANSFERS
    % These are transfers computed on the assumption that Agent 2 cannot
    % borrow and lend. The transfers are the difference between his
    % consumption and after tax earning (l . U_l/U_c)
    Trans=C2-L2.*ul2./uc2;
    
    
     % Income
    AfterTaxWageIncome_Agent2=L2.*ul2./uc2+Trans;
    AfterTaxWageIncome_Agent1=L1.*ul1./uc1+Trans;
    % Gini Coeff
    GiniCoeff=(AfterTaxWageIncome_Agent2 +2*AfterTaxWageIncome_Agent1)./(AfterTaxWageIncome_Agent2+AfterTaxWageIncome_Agent1)-3/2;
   
    

    % UPDATE THE SIMULATION HISTORY
    
    
       RHist(i+1)=RRPrime(sHist(i+1));
    xHist(i+1)=xxPrime(sHist(i+1)) ;
    btildHist(i+1)=btildprime(sHist(i+1)) ;
    TauHist(i+1)=Tau(sHist(i+1));
    YHist(i+1)=y(sHist(i+1));
    TransHist(i+1)=Trans(sHist(i+1));
    c1Hist(i+1)=C1(sHist(i+1));
    c2Hist(i+1)=C2(sHist(i+1));
    l1Hist(i+1)=L1(sHist(i+1));
    l2Hist(i+1)=L2(sHist(i+1));
    IncomeFromAssets_Agent1Hist(i)=-btildHist(i).*(IntHist(i)-1);
    AfterTaxWageIncome_Agent1Hist(i+1)=AfterTaxWageIncome_Agent1(sHist(i+1));
    AfterTaxWageIncome_Agent2Hist(i+1)=AfterTaxWageIncome_Agent2(sHist(i+1));
    Theta_1Hist(i+1)=theta_1(sHist(i+1));
    Theta_2Hist(i+1)=theta_2(sHist(i+1));
    GiniCoeffHist(i+1)=GiniCoeff(sHist(i+1)); 


    
     if mod(i,1000)==0 || i==NumSim-1
        disp('Running Simulation, t=')
        disp(i)
        toc
        tic
SimData.sHist=sHist;
        SimData.Theta_1Hist=Theta_1Hist;
SimData.Theta_2Hist=Theta_2Hist;
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
SimData.GiniCoeffHist=GiniCoeffHist;
save('~/Golosov-Sargent/Data/temp/SimData.mat','SimData')
    end
   
end

end

% BUDGET CONSTRAINTS
% c1Hist(3:5)-btildHist(3:5)-AfterTaxWageIncome_Agent1Hist(3:5)-IncomeFromAssets_Agent1Hist(2:4)
% c2Hist(3:5)-AfterTaxWageIncome_Agent2Hist(3:5)
% YHist(3:5)-n1*c1Hist(3:5)-n2*c2Hist(3:5)-gHist(3:5)
% YHist(3:5)-n1*l1Hist(3:5)*theta_1-n2*l2Hist(3:5)*theta_2
% gHist(3:5)+n1*btildHist(3:5)+TransHist(3:5)-TauHist(3:5).*(n1*theta_1*l1H
% ist(3:5)+n2*theta_2*l2Hist(3:5))-n1*btildHist(2:4).*IntHist(2:4)