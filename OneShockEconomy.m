

SetParaStruc;
 clc
 clear all
 close all
s_=1;
texpath='C:\Users\Anmol\Dropbox\2011RA\FiscalPolicy\GolosovProjectCode\Tom Example\BGP\Tex\Calibration\';






%% Build Grid for the state variables

load('Data/Calibration/cPhMed.mat')
Para.theta_2=0;
   Para.ApproxMethod='spli';
  Para.u2btildGridSize=15;
  Para.RGridSize=15;
  Para.OrderOfAppx_u2btild=10;
  Para.OrderOfApprx_R=10;

Para.u2btildMin=-1;
Para.u2btildMax=3;
Para.RMin=3;
Para.RMax=4;
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
    
    
    figure()
    fplot(@(x) funeval(c0SS(1,:)' ,VSS(1),[x mean(Para.RGrid)]),[Para.u2btildMin Para.u2btildMax])
    figure()
    fplot(@(R) funeval(c0SS(1,:)' ,VSS(1),[mean(Para.u2bdiffGrid) R]),[Para.RMin Para.RMax])
    
    
    
    RGrid.RMin=Para.RMin;
 RGrid.RMax=Para.RMax;
     StationaryPolicyRulesStore=vertcat(squeeze(StationaryPolicyRulesStore(1,:,:)),squeeze(StationaryPolicyRulesStore(2,:,:)));
               
[cHat,VHat]=GetVHat(Para,RGrid,c0SS,VSS,x_state,StationaryPolicyRulesStore);
VHatData=load([Para.datapath '/cVHat.mat'])
Para.StoreFileName='/cVHat.mat'
Para.flagPlot2PeriodDrifts=0
GetPlotsForFinalSolution(Para)

 load([Para.datapath '/cVHat.mat'])
 xState=fsolve(@(x) GetCrossingPoints(x,1,c,V,PolicyRulesStore,x_state,Para),[1 2.5])

 %  funeval(c(1,:)',V(1),[xState])/((GetStationaryValue(xState,Para))/(1-Para.beta))
%  
%  
%  n1=Para.n1;
% n2=Para.n2;
% alpha_1=Para.alpha_1;
% alpha_2=Para.alpha_2;
% theta_1=Para.theta_1;
% theta_2=Para.theta_2;
% psi=Para.psi;
% beta=Para.beta;
% u2btild=xState(1);
% R=xState(2);
%     [PolicyRulesInit]=GetInitialApproxPolicy([u2btild R s_],x_state,PolicyRulesStore);
%     [PolicyRules, V_new,exitflag,~]=CheckGradNAG(u2btild,R,s_,c,V,PolicyRulesInit,Para,0) ;
% 
% %[xState,fvec,xit]=fminunc(@(xState) GetStationaryValue(xState,Para),[1.5 2.8] )
% s=1;
% u2bdiff=xState(1);
% RR=xState(2);
% [x,fvec,exitval]=fsolve(@(x) SteadyStateCharacterization(x,u2bdiff,RR,Para,1) ,[PolicyRules(1:3)]);
% Para.theta=[Para.theta_1 Para.theta_2];
% Para.alpha=[Para.alpha_1 Para.alpha_2];
% Par=Para;
% u2btild=u2bdiff;
% R=RR;
% s_=s;
% n1=Para.n1;
% n2=Para.n2;
% ctol=Para.ctol;
% 
% %% GET THE Policy Rules
% psi= Par.psi;
% beta =  Par.beta;
% P = Par.P;
% theta_1 = Par.theta(1);
% theta_2 = Par.theta(2);
% g = Par.g;
% alpha = Par.alpha;
% sigma = 1;
% c1_1=x(1);
% c1_2=x(2);
% c2_1=x(3);
% 
% %compute components from unconstrained guess
% [c2_2 grad_c2_2] = computeC2_2(c1_1,c1_2,c2_1,R,s_,P,sigma);
% [l1 l1grad l2 l2grad] = computeL(c1_1,c1_2,c2_1,c2_2,grad_c2_2,...
%     theta_1,theta_2,g,n1,n2);
% [btildprime grad_btildprime] = computeBtildeprime(c1_1,c1_2,c2_1,c2_2,grad_c2_2,l1,l2,l1grad,l2grad,...
%     u2btild,s_,psi,beta,P);
% 
% % x' - definition
% u2btildprime=psi*[c2_1^(-1) c2_2^(-1)].*btildprime;
% 
% % State next period
% X(1,:) = [psi*c2_1^(-1)*btildprime(1),c2_1^(-1)/c1_1^(-1)];%state next period
% X(2,:) = [psi*c2_2^(-1)*btildprime(2),c2_2^(-1)/c1_2^(-1)];%state next period
% Vobj = P(s_,1)*(alpha(1)*uBGP(c1_1,l1(1),psi)+alpha(2)*uBGP(c2_1,l2(1),psi));
% Vobj = Vobj + P(s_,2)*(alpha(1)*uBGP(c1_2,l1(2),psi)+alpha(2)*uBGP(c2_2,l2(2),psi));
% 
% PolicyRules=[c1_1 c1_2 c2_1 c2_2 l1(1) l1(2) l2(1) l2(2) btildprime c2_1^(-1)/c1_1^(-1) c2_2^(-1)/c1_2^(-1) u2btildprime(1) u2btildprime(2)];
%     %PolicyRules=[c1_1 c1_2 c2_1 c2_2 l1(1) l1(2) l2(1) l2(2) btildprime c2_1^(-1)/c1_1^(-1) c2_2^(-1)/c1_2^(-1) u2btildprime(1) u2btildprime(2)]
%     c1=PolicyRules(1:2);
%     c2=PolicyRules(3:4);
%     l1=PolicyRules(5:6);
%     l2=PolicyRules(7:8);
%     ul2=(1-psi)./(1-l2);
%     uc2=psi./c2;
%     ul1=(1-psi)./(1-l1);
%     uc1=psi./c1;
%     Rprime=PolicyRules(end-3:end-2);
%     % x' - u_c_2* btildprime
%     u2btildprime=PolicyRules(end-1:end);
%     % btildprime - x'/u_c2
%     btildprime=PolicyRules(9:10);
% 
%     ul2=(1-psi)./(1-l2);
%     uc2=psi./c2;
%     ul1=(1-psi)./(1-l1);
%     uc1=psi./c1;
%     Rprime=PolicyRules(end-3:end-2);
%     % x' - u_c_2* btildprime
%     u2btildprime=PolicyRules(end-1:end);
%     % btildprime - x'/u_c2
%     btildprime=PolicyRules(9:10);
%    
%     % TAU - From the WAGE optimality of Agent 2
%     Tau=1-(ul2./(theta_2.*uc2));
%     
%     % OUTPUT
%     y(1)=c1(1)*n1+c2(1)*n2+g(1);
%     y(2)=c1(2)*n1+c2(2)*n2+g(2);
%     
%     % TRANSFERS
%     % These are transfers computed on the assumption that Agent 2 cannot
%     % borrow and lend. The transfers are the difference between his
%     % consumption and after tax earning (l . U_l/U_c)
%     Trans=c2-l2.*ul2./uc2;
%     
%      % Income
%     AfterTaxWageIncome_Agent2=l2.*ul2./uc2;
%     AfterTaxWageIncome_Agent1=l1.*ul1./uc1;
%     % Gini Coeff
%     GiniCoeff=(AfterTaxWageIncome_Agent2 +2*AfterTaxWageIncome_Agent1)./(AfterTaxWageIncome_Agent2+AfterTaxWageIncome_Agent1)-3/2;
%     Decomp=[Para.g(2)-Para.g(1),Para.n1*(btildprime(2)-btildprime(1)),(Para.n1+Para.n2)*(Trans(2)-Trans(1)),Para.n1*(AfterTaxWageIncome_Agent1(2)-AfterTaxWageIncome_Agent1(1)),Para.n2*(AfterTaxWageIncome_Agent2(2)-AfterTaxWageIncome_Agent2(1))]
%     rowLabels = {'$g\_=g_l$','$g\_=g_h$'};
%     columnLabels = {'$g_h-g_l$','$n_1[b''(h)-b''(l)]$','$(n_1+n_2)[T(h)-T(l)]$','$n_1\theta_1[l_1(h)\tau(h)-l_1(l)\tau(l)]$', '$n_2\theta_2[l_2(h)\tau(h)-l_2(l)\tau(l)]$'};
%     matrix2latex([Decomp;Decomp], [texpath 'GovHighLowPolicyRules.tex'] , 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.4f', 'size', 'tiny');
