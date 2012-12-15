% CAPTION : fig:flagPoints - This plots the sucess of the optimizer to
% solve the FOC at the points selected in the state space for the final set of coeffecients. The red points
% denote failure.
function GetPlotsForFinalSolution(Para,Domain)
flagPlot2PeriodDrifts=Para.flagPlot2PeriodDrifts;
oldplotpath=Para.plotpath;
load([Para.datapath Para.StoreFileName])
close all;
plotpath=oldplotpath;
mkdir(plotpath);
datapath='Data/Calibration/';
disp('Govt Exp')
g=Para.g;
Para.P
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

xSolved=x_state(IndxSolved,:);
xUnSolved=x_state(IndxUnSolved,:);

figure()

scatter(squeeze(xSolved(:,1)),squeeze(xSolved(:,2)),'b','filled')
hold on
scatter(squeeze(xUnSolved(:,1)),squeeze(xUnSolved(:,2)),'r','filled')
hold on
xlabel('$x$','Interpreter','Latex')
ylabel('$R$','Interpreter','Latex')
print(gcf,'-dpng',[plotpath 'flagPoints.png'])

% CAPTION : fig:NumFOCSolved - This plots shows the number of points that
% the FOC had a solution across iterations
figure()
plot(max(cdiff,[],2))
xlabel('Iteration');
ylabel('Max of Coefficient Difference');
print(gcf,'-dpng',[plotpath 'CoeffConvergence.png'])


figure()
plot(ErrorInSupNorm)
xlabel('Iteration');
ylabel('Max of Coefficient Difference');
print(gcf,'-dpng',[plotpath 'CoeffConvergenceSupNorm.png'])

u2btildLL=Para.u2btildMin;
u2btildUL=Para.u2btildMax;
ucbtild_bounds = [u2btildLL,u2btildUL];
Rbounds=[min(Para.RGrid),max(Para.RGrid)];
if nargin==2
ucbtild_bounds=Domain.xBounds;
Rbounds=Domain.RBounds;
end
%Caption : fig:FunctionalConvergence - This figure plots the value function
% with respect to $\tilde{b}_2$ across iterations. The four panels refer to
% vaules of R. The red line is the first iteration

figure()
% Fix s_
s_=1;
RList=linspace(Rbounds(1),Rbounds(2),4);
%for l=1:length(ListIterations)
%    load([ datapath 'c' num2str(ListIterations(l)) '.mat'])

for Rctr=1:4
    subplot(2,2,Rctr)
    fplot(@(u2btild) funeval(c(s_,:)',V(s_),[u2btild RList(Rctr)]),[ucbtild_bounds(1) ucbtild_bounds(2)]);
end
xlabel('$x$','Interpreter','Latex')
title(['$\rho=$' num2str(RList(Rctr))],'Interpreter','Latex')
hold on


print(gcf,'-dpng',[plotpath 'FunctionalConvergence.png'])

%
% % ChebError




numtest=15;
for n=1:numtest
    %
    %
    u2btild=ucbtild_bounds(1)+(ucbtild_bounds(2)-ucbtild_bounds(1))*rand;
    R=Rbounds(1)+(Rbounds(2)-Rbounds(1))*rand;
    
    xTarget(n,:)=[u2btild R s_];
    [PolicyRulesInit]=GetInitialApproxPolicy(xTarget(n,:),x_state,PolicyRulesStore);
    [PolicyRules, V_new,exitflag,~]=CheckGradNAG(u2btild,R,s_,c,V,PolicyRulesInit,Para,0) ;
    VDirect=funeval(c(s_,:)',V(s_),xTarget(n,1:2));
    Check(n)=(VDirect-V_new)/V_new;
    
    %Do optimization
    %Vopt = CheckOpt(u2btild,R,s_,c,V,PolicyRulesInit,Para);
    %Check2(n) = (Vopt-V_new)/V_new;
    %[Vopt1,Vopt2] = CheckOpt(u2btild,R,s_,c,V,PolicyRules(1:3),Para);
    %Check3(n) = (Vopt1-V_new)/V_new;
    %Check4(n) = (Vopt2-V_new)/V_new;
    if ~(exitflag==1)
        colFOC(n,:)=[1 0 0];
    else
        colFOC(n,:)=[0 0 1];
    end
    
    
end

figure()
plot(Check)
xlabel('Number of Test Points')
ylabel('Percentage Error')
print(gcf,'-dpng',[plotpath 'ChebError.png'])

%
% % Caption : fig:ValueFunction - This plot depicts the value function
figure()


figure()

for Rctr=1:4
    subplot(2,2,Rctr)
    fplot(@(u2btild) funeval(c(s_,:)',V(s_),[u2btild RList(Rctr)]),[ucbtild_bounds(1) ucbtild_bounds(2)],'-k');
    xlabel('$x$','Interpreter','Latex')
    title(['$\rho=$' num2str(RList(Rctr))],'Interpreter','Latex')
    hold on
end
print(gcf,'-dpng',[plotpath 'ValueFunctionx.png'])

% % Caption : fig:ValueFunction - This plot depicts the value function
figure()
xlist=linspace(min(Para.u2bdiffGrid),max(Para.u2bdiffGrid),4);
for xctr=1:4
    subplot(2,2,xctr)
    fplot(@(R) funeval(c(s_,:)',V(s_),[xlist(xctr) R]),[Rbounds(1) Rbounds(2)],'-k');
    xlabel('$R$','Interpreter','Latex')
    title(['$x=$' num2str(xlist(xctr))],'Interpreter','Latex')
    hold on
end

print(gcf,'-dpng',[plotpath 'ValueFunctionR.png'])

% % Caption : fig:ValueFunction - This plot depicts the value function
figure()
xlist=0
for xctr=1:1
    fplot(@(R) funeval(c(s_,:)',V(s_),[xlist(xctr) R]),[Rbounds(1) Rbounds(2)],'-k');
    xlabel('$\rho$','Interpreter','Latex')
    title(['$x=$' num2str(xlist(xctr))],'Interpreter','Latex')
    hold on
end

print(gcf,'-dpng',[plotpath 'ValueFunctionRx_0.png'])


%% Policy Rules entire state space
% Caption : fig:PolicyRules - This plot depicts the $\tilde{b}'_2$ as a function of $\tilde{b}_2$
figu2BtildePrime =figure('Name','x');
figEDeltaX=figure('Name','EDeltaX');
figEDeltaR=figure('Name','EDeltaRho');
figVDeltaX=figure('Name','VDeltaX');
figVDeltaR=figure('Name','VDeltaRho');

figBtildePrime =figure('Name','btild');
figRprime=figure('Name','$\rho$');
figFOCRes =figure('Name','FOCRes');
figRR = figure('Name', 'RRGraph');
figHLLH = figure('Name', 'HLLHPlot');
u2bdiffFineGrid=linspace(ucbtild_bounds(1),ucbtild_bounds(2),35);
RFineGrid=linspace(Rbounds(1),Rbounds(2),35);
RList=linspace(Rbounds(1),Rbounds(2),4);
u2bdiffList=linspace(ucbtild_bounds(1),ucbtild_bounds(2),4);
s_=1;
for Rctr=1:4
    for u2btildctr=1:length(u2bdiffFineGrid)
        R=RList(Rctr);
        u2btild=u2bdiffFineGrid(u2btildctr);
        [PolicyRulesInit]=GetInitialApproxPolicy([u2btild R s_] ,x_state,PolicyRulesStore);
        [PolicyRules, V_new,exitflag,fvec]=CheckGradNAG(u2btild,R,s_,c,V,PolicyRulesInit,Para,0);
        if exitflag==1
            IndxPrint(u2btildctr)=1;
        else
            IndxPrint(u2btildctr)=0;
        end
        
        FOCRes(u2btildctr)=max(abs(fvec));
        u2BtildePrime(u2btildctr,:)=PolicyRules(end-1:end);
        BtildePrime(u2btildctr,:)=PolicyRules(end-5:end-4);
        Rprime(u2btildctr,:)=PolicyRules(end-3:end-2)-R;
        EDeltaX(u2btildctr)=sum(Para.P(s_,:).*u2BtildePrime(u2btildctr,:))-u2btild;
        VDeltaX(u2btildctr)=sum(Para.P(s_,:).*(u2BtildePrime(u2btildctr,:)-[u2btild u2btild]).^2)-EDeltaX(u2btildctr).^2;
        EDeltaR(u2btildctr)=sum(Para.P(s_,:).*Rprime(u2btildctr,:))-R;
        VDeltaR(u2btildctr)=sum(Para.P(s_,:).*(Rprime(u2btildctr,:)-[R R]).^2)-EDeltaR(u2btildctr).^2;
        if flagPlot2PeriodDrifts==1
        %Compute u2btildlh
        u2btild = u2BtildePrime(u2btildctr,1);
        R = Rprime(u2btildctr,1)+R;
        [PolicyRulesInit]=GetInitialApproxPolicy([u2btild R 1] ,x_state,PolicyRulesStore);
        [PolicyRules, V_new,exitflag,fvec]=CheckGradNAG(u2btild,R,1,c,V,PolicyRulesInit,Para,0);
        u2btildLHHL(u2btildctr,1) = PolicyRules(end);
        u2btild = u2BtildePrime(u2btildctr,2);
        R = Rprime(u2btildctr,2)+R;
        [PolicyRulesInit]=GetInitialApproxPolicy([u2btild R s_] ,x_state,PolicyRulesStore);
        [PolicyRules, V_new,exitflag,fvec]=CheckGradNAG(u2btild,R,s_,c,V,PolicyRulesInit,Para,0);
        u2btildLHHL(u2btildctr,2) = PolicyRules(end-1);
        end
    end
    
    
    figure(figEDeltaX)
    subplot(2,2,Rctr)
    plot(u2bdiffFineGrid, EDeltaX,'k','LineWidth',2)
    xlabel('$x$','Interpreter','Latex')
    ylabel('EDeltaX','Interpreter','Latex')
    title(['$\rho=$' num2str(RList(Rctr))],'Interpreter','Latex')
    
    
    
    
    
    figure(figVDeltaX)
    subplot(2,2,Rctr)
    plot(u2bdiffFineGrid, VDeltaX,'k','LineWidth',2)
    xlabel('$x$','Interpreter','Latex')
    ylabel('VDeltaX','Interpreter','Latex')
    title(['$\rho=$' num2str(RList(Rctr))],'Interpreter','Latex')
    
    
    
    figure(figEDeltaR)
    subplot(2,2,Rctr)
    plot(u2bdiffFineGrid, EDeltaR,'k','LineWidth',2)
    xlabel('$x$','Interpreter','Latex')
    ylabel('EDeltaR','Interpreter','Latex')
    title(['$\rho=$' num2str(RList(Rctr))],'Interpreter','Latex')
    
    
    
    figure(figVDeltaR)
    subplot(2,2,Rctr)
    plot(u2bdiffFineGrid, VDeltaR,'k','LineWidth',2)
    xlabel('$x$','Interpreter','Latex')
    ylabel('VDeltaR','Interpreter','Latex')
    title(['$\rho=$' num2str(RList(Rctr))],'Interpreter','Latex')
    
    
    
    figure(figFOCRes)
    subplot(2,2,Rctr)
    plot(u2bdiffFineGrid, FOCRes,'k','LineWidth',2)
    if Rctr==1
        legend('g_l','g_h')
    end
    xlabel('$x$','Interpreter','Latex')
    ylabel('FOCRes','Interpreter','Latex')
    title(['$\rho=$' num2str(RList(Rctr))],'Interpreter','Latex')
    
    figure(figBtildePrime)
    subplot(2,2,Rctr)
    plot(u2bdiffFineGrid(logical(IndxPrint)), BtildePrime(logical(IndxPrint),1),'k','LineWidth',2)
    hold on
    plot(u2bdiffFineGrid(logical(IndxPrint)), BtildePrime(logical(IndxPrint),2),':k','LineWidth',2)
    hold on
    if Rctr==1
        legend('g_l','g_h')
    end
    %plot(u2bdiffFineGrid, 0*u2bdiffFineGrid,':k','LineWidth',2);
    hold on
    %plot(u2bdiffFineGrid,repmat([u2btildLL u2btildUL],length(u2bdiffFineGrid),1)-[u2bdiffFineGrid' u2bdiffFineGrid'] ,':r')
    %
    xlabel('$x$','Interpreter','Latex')
    ylabel('$\tilde{b}_2$','Interpreter','Latex')
    title(['$\rho=$' num2str(RList(Rctr))],'Interpreter','Latex')
    
    figure(figu2BtildePrime)
    subplot(2,2,Rctr)
    plot(u2bdiffFineGrid(logical(IndxPrint)), u2BtildePrime(logical(IndxPrint),1)- u2bdiffFineGrid(logical(IndxPrint))','k','LineWidth',2)
    hold on
    plot(u2bdiffFineGrid(logical(IndxPrint)), u2BtildePrime(logical(IndxPrint),2)-u2bdiffFineGrid(logical(IndxPrint))',':k','LineWidth',2)
    hold on
    if Rctr==1
        legend('g_l','g_h')
    end
    %plot(u2bdiffFineGrid, 0*u2bdiffFineGrid,':k','LineWidth',2);
    hold on
    %plot(u2bdiffFineGrid,repmat([u2btildLL u2btildUL],length(u2bdiffFineGrid),1)-[u2bdiffFineGrid' u2bdiffFineGrid'] ,':r')
    %
    xlabel('$x$','Interpreter','Latex')
    ylabel('$x(s)-x$','Interpreter','Latex')
    title(['$\rho=$' num2str(RList(Rctr))],'Interpreter','Latex')
    if flagPlot2PeriodDrifts==1
    figure(figHLLH)
    subplot(2,2,Rctr)
    plot(u2bdiffFineGrid, u2btildLHHL(:,1)'-u2bdiffFineGrid,'k','LineWidth',2);
    hold on
    plot(u2bdiffFineGrid, u2btildLHHL(:,2)'-u2bdiffFineGrid,':k','LineWidth',2);
    hold on
    if Rctr==1
        legend('LH','HL')
    end
    xlabel('$x$','Interpreter','Latex')
    ylabel('$x_{t+2}-x_t$','Interpreter','Latex')
    title(['$\rho=$' num2str(RList(Rctr))],'Interpreter','Latex')
    end
    
    %
    figure(figRprime)
    subplot(2,2,Rctr)
    plot(u2bdiffFineGrid(logical(IndxPrint)), Rprime(logical(IndxPrint),1),'k','LineWidth',2);
    hold on
    plot(u2bdiffFineGrid(logical(IndxPrint)), Rprime(logical(IndxPrint),2),':k','LineWidth',2);
    xlabel('$x$','Interpreter','Latex')
    ylabel('$\rho(s)-\rho$','Interpreter','Latex')
    title(['$\rho=$' num2str(RList(Rctr))],'Interpreter','Latex')
end


for u2btildctr=1:4
    for Rctr=1:length(RFineGrid)
        R=RFineGrid(Rctr);
        u2btild=u2bdiffList(u2btildctr);
        [PolicyRulesInit]=GetInitialApproxPolicy([u2btild R s_] ,x_state,PolicyRulesStore);
        [PolicyRules, V_new,exitflag,fvec]=CheckGradNAG(u2btild,R,s_,c,V,PolicyRulesInit,Para,0);
        if exitflag==1
            IndxPrint(u2btildctr)=1;
        else
            IndxPrint(u2btildctr)=0;
        end
        
        FOCRes(Rctr)=max(abs(fvec));
        u2BtildePrime(Rctr,:)=PolicyRules(end-1:end);
        BtildePrime(Rctr,:)=PolicyRules(end-5:end-4);
        Rprime(Rctr,:)=PolicyRules(end-3:end-2)-R;
        EDeltaX(Rctr)=sum(Para.P(s_,:).*u2BtildePrime(Rctr,:))-u2btild;
        VDeltaX(Rctr)=sum(Para.P(s_,:).*(u2BtildePrime(Rctr,:)-[u2btild u2btild]).^2)-EDeltaX(Rctr).^2;
        EDeltaR(Rctr)=sum(Para.P(s_,:).*Rprime(Rctr,:))-R;
        VDeltaR(Rctr)=sum(Para.P(s_,:).*(Rprime(Rctr,:)-[R R]).^2)-EDeltaR(Rctr).^2;
        
    end
    figure(figRR)
    subplot(2,2,u2btildctr);
    hold on
    plot(RFineGrid,Rprime(:,1)','k','LineWidth',2);
    hold on
    plot(RFineGrid,Rprime(:,2)',':k','LineWidth',2);
    xlabel('$\rho$','Interpreter','Latex')
    ylabel('$\rho(s)-\rho$','Interpreter','Latex')
    title(['$x=$' num2str(u2bdiffList(u2btildctr))],'Interpreter','Latex')
end

print(figFOCRes,'-dpng',[plotpath 'FOCResFullDomain.png'])

print(figEDeltaX,'-dpng',[plotpath 'EDeltaX.png'])
print(figEDeltaR,'-dpng',[plotpath 'EDeltaR.png'])
print(figVDeltaX,'-dpng',[plotpath 'VDeltaX.png'])
print(figVDeltaR,'-dpng',[plotpath 'VDeltaR.png'])

print(figu2BtildePrime,'-dpng',[plotpath 'u2BtildePrimeFullDomain.png'])
print(figBtildePrime,'-dpng',[plotpath 'BtildePrimeFullDomain.png'])
print(figRR,'-dpng',[plotpath 'RPrimeFullDomain.png'])


% fixed point policy rules
figDeltaXRSS=figure('Name','Steady State Policyies');
[ RSS,xSS,~ ] = findSteadyState( 0,mean(Para.RGrid),Para);
% Change in x
for u2btildctr=1:length(u2bdiffFineGrid)
        R=RSS;
        u2btild=u2bdiffFineGrid(u2btildctr);
        [PolicyRulesInit]=GetInitialApproxPolicy([u2btild R s_] ,x_state,PolicyRulesStore);
        [PolicyRules, V_new,exitflag,fvec]=CheckGradNAG(u2btild,R,s_,c,V,PolicyRulesInit,Para,0);
        FOCRes(u2btildctr)=max(abs(fvec));
        u2BtildePrime(u2btildctr,:)=PolicyRules(end-1:end);
        BtildePrime(u2btildctr,:)=PolicyRules(end-5:end-4);
        Rprime(u2btildctr,:)=PolicyRules(end-3:end-2)-R;
        EDeltaX(u2btildctr)=sum(Para.P(s_,:).*u2BtildePrime(u2btildctr,:))-u2btild;
        VDeltaX(u2btildctr)=sum(Para.P(s_,:).*(u2BtildePrime(u2btildctr,:)-[u2btild u2btild]).^2)-EDeltaX(u2btildctr).^2;
        EDeltaR(u2btildctr)=sum(Para.P(s_,:).*Rprime(u2btildctr,:))-R;
        VDeltaR(u2btildctr)=sum(Para.P(s_,:).*(Rprime(u2btildctr,:)-[R R]).^2)-EDeltaR(u2btildctr).^2;
    
end
 figure(figDeltaXRSS)
   subplot(1,2,1)
    plot(u2bdiffFineGrid(logical(IndxPrint)), u2BtildePrime(logical(IndxPrint),1)- u2bdiffFineGrid(logical(IndxPrint))','k','LineWidth',2)
    hold on
    plot(u2bdiffFineGrid(logical(IndxPrint)), u2BtildePrime(logical(IndxPrint),2)-u2bdiffFineGrid(logical(IndxPrint))',':k','LineWidth',2)
    hold on
        legend('g_l','g_h')
  xlabel('x')
    ylabel('x(s)-x')
    title(['\rho=' RSS])
    

    for Rctr=1:length(RFineGrid)
        R=RFineGrid(Rctr);
        u2btild=xSS;
        [PolicyRulesInit]=GetInitialApproxPolicy([u2btild R s_] ,x_state,PolicyRulesStore);
        [PolicyRules, V_new,exitflag,fvec]=CheckGradNAG(u2btild,R,s_,c,V,PolicyRulesInit,Para,0);
    
        FOCRes(Rctr)=max(abs(fvec));
        u2BtildePrime(Rctr,:)=PolicyRules(end-1:end);
        BtildePrime(Rctr,:)=PolicyRules(end-5:end-4);
        Rprime(Rctr,:)=PolicyRules(end-3:end-2)-R;
        EDeltaX(Rctr)=sum(Para.P(s_,:).*u2BtildePrime(Rctr,:))-u2btild;
        VDeltaX(Rctr)=sum(Para.P(s_,:).*(u2BtildePrime(Rctr,:)-[u2btild u2btild]).^2)-EDeltaX(Rctr).^2;
        EDeltaR(Rctr)=sum(Para.P(s_,:).*Rprime(Rctr,:))-R;
        VDeltaR(Rctr)=sum(Para.P(s_,:).*(Rprime(Rctr,:)-[R R]).^2)-EDeltaR(Rctr).^2;
        
    end
    figure(figDeltaXRSS)
    subplot(1,2,2);
    plot(RFineGrid,Rprime(:,1)','k','LineWidth',2);
    hold on
    plot(RFineGrid,Rprime(:,2)',':k','LineWidth',2);
    xlabel('\rho')
    ylabel('\rho(s)-\rho')
    title(['x=' xSS])

       legend('g_l','g_h')
print(figDeltaXRSS,'-dpng',[plotpath 'figDeltaXRSS.png'])
    