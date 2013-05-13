function  [SimData]=RunSimulationsFromT1AltThetaShocksUsingPolicyRules(Coeff_xhat,Coeff_Rhat,xhat,Rhat,x0,R0,NumSim,Para,rHist0)
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

 disp('theta Exp')
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


% RUN SIMULATION
Theta_1Hist=zeros(NumSim,1);
Theta_2Hist=zeros(NumSim,1);
xHist=zeros(NumSim,1);
RHist=zeros(NumSim,1);
sHist(1)=1;


% INITIALIZE - t=0
xHist(1)=x0;
RHist(1)=R0;
tic

for i=1:NumSim-1
    % ------STATE (t) - x,R,s_ ------------------------------------------
    x=xHist(i);
    R=RHist(i);
    
    % DRAW THE s' ~ P(s,:) if flagUseExistingShocks is set to no
    if strcmpi(flagUseExistingShocks,'yes')
        sHist(i+1)=sum(~(CumP(sHist(i),:)-rHist0(i+1)>0));

    else
            sHist(i+1)=sum(~(CumP(sHist(i),:)-rand>0));
    end
    % UPDATE THE SIMULATION HISTORY
    
    RHist(i+1)=funeval(Coeff_Rhat(sHist(i+1),:)',Rhat(sHist(i+1)),[x R]);
    xHist(i+1)=funeval(Coeff_xhat(sHist(i+1),:)',xhat(sHist(i+1)),[x R]);
    Theta_1Hist(i+1)=theta_1(sHist(i+1));
    Theta_2Hist(i+1)=theta_2(sHist(i+1));

    
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