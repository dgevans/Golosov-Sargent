% -------------------------------------------------------------------------
% This script uses the policy rules to create a linear and non linear approximation of LOM
% of state variables  - x,R
%--------------------------------------------------------------------------

% Non Linear : x(s)-xSS=B(x-xSS,R-RSS,s) and  R(s)-RSS=B(x-xSS,R-RSS,s)
% Linear
%{ 

| x(s) - xSS |               | x-xSS    |
|            |    =  B(s)    |          |
| R(s) - RSS |               | R -RSS   | 


%}



clear all
clc
close all
% Load the coeff
load('Data/temp/csigmaMed.mat')
% Find the SS
% Using the crossing point of policy rules
 xState=fsolve(@(x) GetCrossingPoints(x,1,c,V,PolicyRulesStore,x_state,Para),[0 mean(Para.RGrid)]);
x0=xState(1);
R0=xState(2);
% Checking against david's code
[ x,R,PolicyRule ] = findSteadyState( x0,R0,Para);
% Approximation 1 : Non Linear Projection

% Setup the approximation domain
xhatMin=-(x0-min(Para.u2bdiffGrid))/2;
xhatMax=(max(Para.u2bdiffGrid)-x0)/2;
RhatMin=-(R0-min(Para.RGrid))/3;
RhatMax=(max(Para.RGrid)-R0)/3;
xhatGridSize=length(Para.u2bdiffGrid);
RhatGridSize=length(Para.u2bdiffGrid);
xhat(1) = fundefn(Para.ApproxMethod,[Para.OrderOfAppx_u2btild Para.OrderOfApprx_R ] ,[xhatMin RhatMin],[xhatMax RhatMax]);
xhat(2) = fundefn(Para.ApproxMethod,[Para.OrderOfAppx_u2btild Para.OrderOfApprx_R ] ,[xhatMin RhatMin],[xhatMax RhatMax]);

Rhat(1) = fundefn(Para.ApproxMethod,[Para.OrderOfAppx_u2btild Para.OrderOfApprx_R ] ,[xhatMin RhatMin],[xhatMax RhatMax]);
Rhat (2)= fundefn(Para.ApproxMethod,[Para.OrderOfAppx_u2btild Para.OrderOfApprx_R ] ,[xhatMin RhatMin],[xhatMax RhatMax]);

xhatGrid=linspace(xhatMin,xhatMax,xhatGridSize);
RhatGrid=linspace(RhatMin,RhatMax,RhatGridSize);
s_=1;

% Compute the fitted poits
ctr=1;
for xind=1:xhatGridSize
    for Rind=1:RhatGridSize
                R=R0+RhatGrid(Rind);
        u2btild=x0+xhatGrid(xind);
        ApproximationDomain(ctr,:)=[xhatGrid(xind) RhatGrid(Rind)];
        [PolicyRulesInit]=GetInitialApproxPolicy([u2btild R s_] ,x_state,PolicyRulesStore);
        [PolicyRules, V_new,exitflag,fvec]=CheckGradNAG(u2btild,R,s_,c,V,PolicyRulesInit,Para,0);
        if exitflag==1
            IndxPrint(ctr)=1;
        else
            IndxPrint(ctr)=0;
        end
         Rhatprime(ctr,:)=PolicyRules(end-3:end-2)-[R0 R0];
    xhatprime(ctr,:)=PolicyRules(end-1:end)-[x0 x0];
    
        
ctr=ctr+1;
    end
end


for s=1:2
        Coeff_xhat(s,:)=funfitxy(xhat(s),ApproximationDomain, xhatprime(:,s));
        Coeff_Rhat(s,:)=funfitxy(Rhat(s),ApproximationDomain, Rhatprime(:,s));
        B(s).Val(1,1)=funeval(Coeff_xhat(s,:)',xhat(s),[0 0] ,[1 0]);
        B(s).Val(1,2)=funeval(Coeff_xhat(s,:)',xhat(s),[0 0] ,[0 1]);
        B(s).Val(2,1)=funeval(Coeff_Rhat(s,:)',Rhat(s),[0 0] ,[1 0]);
        B(s).Val(2,2)=funeval(Coeff_Rhat(s,:)',Rhat(s),[0 0] ,[0 1]);
       ApproxError(s).Val= (B(s).Val*ApproximationDomain')'-[xhatprime(:,s) Rhatprime(:,s)];
end



% Store the B(s) matrices
save('Data/temp/PolicyRulesApproximation.mat','xhat','Coeff_xhat','Rhat','Coeff_Rhat','B','ApproxError')
% Check the approximation errors
RhatList=linspace(RhatMin,RhatMax,4);
figure()
for l=1:4
subplot(2,2,l)
[~,FindClosestR]=min(abs(ApproximationDomain(:,2)-RhatList(l)));
xPlotGridIndx=find(ApproximationDomain(:,2)==ApproximationDomain(FindClosestR,2));
plot(ApproximationDomain(xPlotGridIndx,1),xhatprime(xPlotGridIndx,1)-[ApproximationDomain(xPlotGridIndx,1)],'k','LineWidth',2)
hold on
plot(ApproximationDomain(xPlotGridIndx,1),xhatprime(xPlotGridIndx,1)+ApproxError(1).Val(xPlotGridIndx,1)-[ApproximationDomain(xPlotGridIndx,1)],'r','LineWidth',1)
hold on
plot(ApproximationDomain(xPlotGridIndx,1),xhatprime(xPlotGridIndx,2)-[ApproximationDomain(xPlotGridIndx,1)],':k','LineWidth',2)
hold on
plot(ApproximationDomain(xPlotGridIndx,1),xhatprime(xPlotGridIndx,2)+ApproxError(2).Val(xPlotGridIndx,1)-[ApproximationDomain(xPlotGridIndx,1)],':r','LineWidth',1)
title(['R-RSS=' num2str(ApproximationDomain(FindClosestR,2))])
xlabel('x-x_{SS}')
ylabel('x(s)-x')
end

xhatList=linspace(xhatMin,xhatMax,4);
figure()
for l=1:4
subplot(2,2,l)
[~,FindClosestx]=min(abs(ApproximationDomain(:,1)-xhatList(l)));
RPlotGridIndx=find(ApproximationDomain(:,1)==ApproximationDomain(FindClosestx,1));
plot(ApproximationDomain(RPlotGridIndx,2),Rhatprime(RPlotGridIndx,1)-[ApproximationDomain(RPlotGridIndx,2)],'k','LineWidth',2)
hold on
plot(ApproximationDomain(RPlotGridIndx,2),Rhatprime(RPlotGridIndx,1)+ApproxError(1).Val(RPlotGridIndx,2)-[ApproximationDomain(RPlotGridIndx,2)],'r','LineWidth',1)
hold on
plot(ApproximationDomain(RPlotGridIndx,2),Rhatprime(RPlotGridIndx,2)-[ApproximationDomain(RPlotGridIndx,2)],':k','LineWidth',2)
hold on
plot(ApproximationDomain(RPlotGridIndx,2),Rhatprime(RPlotGridIndx,2)+ApproxError(2).Val(xPlotGridIndx,2)-[ApproximationDomain(RPlotGridIndx,2)],':r','LineWidth',1)
title(['x-xSS=' num2str(ApproximationDomain(FindClosestx,1))])
xlabel('R-R_{SS}')
ylabel('R(s)-R')
end