
% % Set the Parallel Config
% err=[];
% try
%     matlabpool('size')
% catch err
% end
% if isempty(err)
%     
%     
%     if(matlabpool('size') == 0)
%         matlabpool close
%            matlabpool open local;
%  
%     end
%     
%     
% end

SetParaStruc;
 clc
 clear all
 close all
s_=1;
load('Data/Calibration/cWar.mat')
texpath='C:\Users\Anmol\Dropbox\2011RA\FiscalPolicy\GolosovProjectCode\Tom Example\BGP\Tex\Calibration\';






%% Build Grid for the state variables
u2btildMin=-3;
u2btildMax=3;
u2btildGrid=linspace(u2btildMin,u2btildMax,Para.u2btildGridSize);

Para.u2bdiffGrid=u2btildGrid;
Para.u2btildLL=u2btildMin;
Para.u2btildUL=u2btildMax;  

RMin=Para.RMin;
RMax=Para.RMax;
RMin=2.2;
RMax=3;
RGrid=linspace(RMin,RMax,Para.RGridSize);
Para.RGrid=RGrid;
GridSize=Para.u2btildGridSize*Para.RGridSize*Para.sSize;
Para.GridSize=GridSize;
Para.u2btildMin=u2btildMin;
Para.u2btildMax=u2btildMax;
Para.RMax=RMax;
Para.RMin=RMin;
%% Define the funtional space
VSS(1) = fundefn(Para.ApproxMethod,[Para.OrderOfAppx_u2btild Para.OrderOfApprx_R ] ,[u2btildMin RMin],[u2btildMax RMax]);
VSS(2) = VSS(1);
GridPoints=[Para.u2btildLL Para.u2btildUL;RMin RMax];
    
    %% INITIALIZE THE COEFF
    %  This function computes c1,c2,l1,l2 and the value for an arbitrary x, R.
    % This section solves for V i.e the value function at the end of period 1
    % with g_t=g for all t >1. since the value function is static we need to
    % solve a equation in c_1 for each x,R. Th function getValueC1 does the job
    % by solving for the two roots of this equation and using the one that
    % supports the highest utility
    tic
    for s_=1:Para.sSize
        n=1;
        if s_==1
            
            
            for u2btildctr=1:Para.u2btildGridSize
                for Rctr=1:Para.RGridSize
                    u2btild_=u2btildGrid(u2btildctr);
                    R_=RGrid(Rctr);
                                         x_state_(n,:)=[u2btild_ R_];

                    %if R_>Rbar(u2btildctr)
                    V0(s_,n)=(GetStationaryValue([u2btild_,R_],Para))/(1-Para.beta);                  
                    n=n+1;
                    %end
                    
                end
            end
             else
            V0(s_,:)=V0(s_-1,:);
            
        end
    end
    
c0SS(1,:)=funfitxy(VSS(1),x_state_,V0(1,:)' );
c0SS(2,:)=c0SS(1,:);
      tic
           [res]=StationaryFOC([1.5 2.808],1,c0SS,VSS,Para)
       
[xx,fvec,exitflag]=fsolve(@(x) StationaryFOC(x,1,c0SS,VSS,Para), [1.5 2.808])
 
%load('Data/Calibration/cWar.mat')
 xState=fsolve(@(x) GetCrossingPoints(x,1,c,V,PolicyRulesStore,x_state,Para),[1.5 2.8])
 
 funeval(c(1,:)',V(1),[xState])/((GetStationaryValue(xState,Para))/(1-Para.beta))
 n1=Para.n1;
n2=Para.n2;
alpha_1=Para.alpha_1;
alpha_2=Para.alpha_2;
theta_1=Para.theta_1;
theta_2=Para.theta_2;
psi=Para.psi;
beta=Para.beta;
u2btild=xState(1);
R=xState(2);
    [PolicyRulesInit]=GetInitialApproxPolicy([u2btild R s_],x_state,PolicyRulesStore);
    [PolicyRules, V_new,exitflag,~]=CheckGradNAG(u2btild,R,s_,c,V,PolicyRulesInit,Para,0) ;

%[xState,fvec,xit]=fminunc(@(xState) GetStationaryValue(xState,Para),[1.5 2.8] )
s=1;
u2bdiff=xState(1);
RR=xState(2);
[x,fvec,exitval]=fsolve(@(x) SteadyStateCharacterization(x,u2bdiff,RR,Para,1) ,[PolicyRules(1:3)]);
Para.theta=[Para.theta_1 Para.theta_2];
Para.alpha=[Para.alpha_1 Para.alpha_2];
Par=Para;
u2btild=u2bdiff;
R=RR;
s_=s;
n1=Para.n1;
n2=Para.n2;
ctol=Para.ctol;

%% GET THE Policy Rules
psi= Par.psi;
beta =  Par.beta;
P = Par.P;
theta_1 = Par.theta(1);
theta_2 = Par.theta(2);
g = Par.g;
alpha = Par.alpha;
sigma = 1;
c1_1=x(1);
c1_2=x(2);
c2_1=x(3);

%compute components from unconstrained guess
[c2_2 grad_c2_2] = computeC2_2(c1_1,c1_2,c2_1,R,s_,P,sigma);
[l1 l1grad l2 l2grad] = computeL(c1_1,c1_2,c2_1,c2_2,grad_c2_2,...
    theta_1,theta_2,g,n1,n2);
[btildprime grad_btildprime] = computeBtildeprime(c1_1,c1_2,c2_1,c2_2,grad_c2_2,l1,l2,l1grad,l2grad,...
    u2btild,s_,psi,beta,P);

% x' - definition
u2btildprime=psi*[c2_1^(-1) c2_2^(-1)].*btildprime;

% State next period
X(1,:) = [psi*c2_1^(-1)*btildprime(1),c2_1^(-1)/c1_1^(-1)];%state next period
X(2,:) = [psi*c2_2^(-1)*btildprime(2),c2_2^(-1)/c1_2^(-1)];%state next period
Vobj = P(s_,1)*(alpha(1)*uBGP(c1_1,l1(1),psi)+alpha(2)*uBGP(c2_1,l2(1),psi));
Vobj = Vobj + P(s_,2)*(alpha(1)*uBGP(c1_2,l1(2),psi)+alpha(2)*uBGP(c2_2,l2(2),psi));

PolicyRules=[c1_1 c1_2 c2_1 c2_2 l1(1) l1(2) l2(1) l2(2) btildprime c2_1^(-1)/c1_1^(-1) c2_2^(-1)/c1_2^(-1) u2btildprime(1) u2btildprime(2)];
    %PolicyRules=[c1_1 c1_2 c2_1 c2_2 l1(1) l1(2) l2(1) l2(2) btildprime c2_1^(-1)/c1_1^(-1) c2_2^(-1)/c1_2^(-1) u2btildprime(1) u2btildprime(2)]
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

    ul2=(1-psi)./(1-l2);
    uc2=psi./c2;
    ul1=(1-psi)./(1-l1);
    uc1=psi./c1;
    Rprime=PolicyRules(end-3:end-2);
    % x' - u_c_2* btildprime
    u2btildprime=PolicyRules(end-1:end);
    % btildprime - x'/u_c2
    btildprime=PolicyRules(9:10);
   
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
    
     % Income
    AfterTaxWageIncome_Agent2=l2.*ul2./uc2;
    AfterTaxWageIncome_Agent1=l1.*ul1./uc1;
    % Gini Coeff
    GiniCoeff=(AfterTaxWageIncome_Agent2 +2*AfterTaxWageIncome_Agent1)./(AfterTaxWageIncome_Agent2+AfterTaxWageIncome_Agent1)-3/2;
    Decomp=[Para.g(2)-Para.g(1),Para.n1*(btildprime(2)-btildprime(1)),(Para.n1+Para.n2)*(Trans(2)-Trans(1)),Para.n1*(AfterTaxWageIncome_Agent1(2)-AfterTaxWageIncome_Agent1(1)),Para.n2*(AfterTaxWageIncome_Agent2(2)-AfterTaxWageIncome_Agent2(1))]
    rowLabels = {'$g\_=g_l$','$g\_=g_h$'};
    columnLabels = {'$g_h-g_l$','$n_1[b''(h)-b''(l)]$','$(n_1+n_2)[T(h)-T(l)]$','$n_1\theta_1[l_1(h)\tau(h)-l_1(l)\tau(l)]$', '$n_2\theta_2[l_2(h)\tau(h)-l_2(l)\tau(l)]$'};
    matrix2latex([Decomp;Decomp], [texpath 'GovHighLowPolicyRules.tex'] , 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.4f', 'size', 'tiny');
