function [ PolicyRules,Para,xState] = SolveSSForOneAgentOneShockEconomy(theta_1,alpha_1,u2btildMin,u2btildMax,RMin,RMax)


%% Build Grid for the state variables

 clc
 close all

 
% Set up the para structure  
SetParaStruc
theta_1old=3.3; % theta high
theta_2old=1;  % theta low
g_l_y=.12; % g low
g_h_y=.13; % g high
n1=1;
n2=1;
tau=.2;
g_Y=mean([g_l_y g_h_y]);
AvfFETarget=.5;
x=fsolve(@(x) GetCalibrationFrischElasticity (x,AvfFETarget,theta_1old,theta_2old,tau,g_Y,n1,n2), [1 1 ]);
gamma=x(1);
Y=x(2);
g=g_Y*Y;
psi=1/(1+gamma);
beta=.9;
%theta_2=0;  % theta low
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
 

Para.theta_2=theta_2;
Para.theta_1=theta_1;
alpha_2=1-alpha_1;
Para.alpha_1=alpha_1*Para.n1;
Para.alpha_2=alpha_2*Para.n2;

  Para.ApproxMethod='spli';
  Para.u2btildGridSize=30;
  Para.RGridSize=30;
  Para.OrderOfAppx_u2btild=29;
  Para.OrderOfApprx_R=29;

Para.u2btildMin=u2btildMin;
Para.u2btildMax=u2btildMax;
Para.RMin=RMin;
Para.RMax=RMax;
Para.u2bdiffGrid=linspace(Para.u2btildMin,Para.u2btildMax,Para.u2btildGridSize);
%Para.u2bdiffGrid=linspace(Para.u2btildMax,Para.u2btildMin,Para.u2btildGridSize);
Para.RGrid=linspace(Para.RMin,Para.RMax,Para.RGridSize);
%% Define the funtional space
VSS(1) = fundefn(Para.ApproxMethod,[Para.OrderOfAppx_u2btild Para.OrderOfApprx_R ] ,[Para.u2btildMin Para.RMin],[Para.u2btildMax Para.RMax]);
VSS(2) = VSS(1);
    u2btildGrid=Para.u2bdiffGrid;
    RGrid=Para.RGrid;
    %% INITIALIZE THE COEFF
    %  This function computes c1,c2,l1,l2 and the value for an arbitrary x, R.
    % This section solves for V i.e the value function at the end of period 1
    % with g_t=g for all t >1. since the value function is static we need to
    % solve a equation in c_1 for each x,R. Th function getValueC1 does the job
    % by solving for the two roots of this equation and using the one that
    % supports the highest utility
    tic
    xInit=[1 1 mean(Para.RGrid)^-1];
    for s_=1:Para.sSize
        n=1;
            
                        for u2btildctr=1:Para.u2btildGridSize
                for Rctr=1:Para.RGridSize
                    u2btild_=u2btildGrid(u2btildctr);
                    R_=RGrid(Rctr);
                                         x_state_(n,:)=[u2btild_ R_];

                    %if R_>Rbar(u2btildctr)
                    res=SolveNoShockProblem([u2btild_,R_ s_],Para);
                    
                    V0(s_,n)=res.Value;
                    StationaryPolicyRulesStore(s_,n,:)=res.PolicyRules;
                    ExitFlag(s_,n)=res.exitflag;
                    if ~(res.exitflag==1)
                        disp([u2btild_ R_])
                    end
                       xInit=res.PolicyRules;
                      xInit=res.PolicyRules;
                    n=n+1;
                    %end
                    
                end
            end
             
        end
    
    
   disp('points where the stationary policies could not be found')
x_state_(logical(~(ExitFlag(1,:)==1)),:)

c0SS(1,:)=funfitxy(VSS(1),x_state_(logical(ExitFlag(1,:)==1),:),V0(1,logical(ExitFlag(1,:)==1))' );
c0SS(2,:)=funfitxy(VSS(2),x_state_(logical(ExitFlag(2,:)==1),:),V0(2,logical(ExitFlag(2,:)==1))' );

    x_state=vertcat([(x_state_) ones(length(x_state_),1)] ,[(x_state_) 2*ones(length(x_state_),1)]);
    
    
    
    RGrid=[];
    RGrid.RMin=RMin;
 RGrid.RMax=RMax;
     StationaryPolicyRulesStore=vertcat(squeeze(StationaryPolicyRulesStore(1,:,:)),squeeze(StationaryPolicyRulesStore(2,:,:)));
               
[cHat,VHat]=GetVHat(Para,RGrid,c0SS,VSS,x_state,StationaryPolicyRulesStore);
%Para.StoreFileName='/cVHat.mat'
%Para.flagPlot2PeriodDrifts=0
%GetPlotsForFinalSolution(Para)

 load([Para.datapath '/cVHat.mat'])
 %options=optimset('TolX',1e-7,'TolFun',1e-7);
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
    % These are transfers computed on the assumption that Agent 2 cannot
    % borrow and lend. The transfers are the difference between his
    % consumption and after tax earning (l . U_l/U_c)
    Trans=c2;
    
    
     % Income
    AfterTaxWageIncome_Agent2=Trans;
    AfterTaxWageIncome_Agent1=l1.*ul1./uc1+Trans;
    % Gini Coeff
    GiniCoeff=(AfterTaxWageIncome_Agent2 +2*AfterTaxWageIncome_Agent1)./(AfterTaxWageIncome_Agent2+AfterTaxWageIncome_Agent1)-3/2;

end

