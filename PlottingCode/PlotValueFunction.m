% CAPTION : fig:flagPoints - This plots the sucess of the optimizer to
% solve the FOC at the points selected in the state space for the final set of coeffecients. The red points
% denote failure.
function PlotValueFunction(Para,Domain)
oldplotpath=Para.plotpath;
load([Para.datapath Para.StoreFileName])
close all;
plotpath=oldplotpath;
mkdir(plotpath);
datapath='Data/Calibration/';
g=Para.g;
S={'-k',':k','.-k'};
Para.P
Para.P;
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
S=length(Para.g);
xSolved=domain(IndxSolved,:);
xUnSolved=domain(IndxUnSolved,:);
gShockNames={'g_l','g_m','g_h'};


figure()

scatter(squeeze(xSolved(:,1)),squeeze(xSolved(:,2)),'b','filled')
hold on
scatter(squeeze(xUnSolved(:,1)),squeeze(xUnSolved(:,2)),'r','filled')
hold on
xlabel('$x$','Interpreter','Latex')
ylabel('$R$','Interpreter','Latex')

% CAPTION : fig:NumFOCSolved - This plots shows the number of points that
% the FOC had a solution across iterations
figure()
plot(max(cdiff,[],2))
xlabel('Iteration');
ylabel('Max of Coefficient Difference');


figure()
plot(ErrorInSupNorm)
xlabel('Iteration');
ylabel('Max of Coefficient Difference');

xLL=Para.xMin;
xUL=Para.xMax;
ucbtild_bounds = [xLL,xUL];
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
    fplot(@(x) funeval(c(s_,:)',V(s_),[x RList(Rctr)]),[ucbtild_bounds(1) ucbtild_bounds(2)]);
end
xlabel('$x$','Interpreter','Latex')
title(['$\rho=$' num2str(RList(Rctr))],'Interpreter','Latex')
hold on


%
% % ChebError




numtest=15;
for n=1:numtest
    %
    %
    x=ucbtild_bounds(1)+(ucbtild_bounds(2)-ucbtild_bounds(1))*rand;
    R=Rbounds(1)+(Rbounds(2)-Rbounds(1))*rand;
    
    xTarget(n,:)=[x R s_];
    [PolicyRulesInit]=GetInitialApproxPolicy(xTarget(n,:),domain,PolicyRulesStore);
    [PolicyRules, V_new,exitflag,~]=CheckGradNAG(x,R,s_,c,V,PolicyRulesInit,Para) ;
    VDirect=funeval(c(s_,:)',V(s_),xTarget(n,1:2));
    Check(n)=(VDirect-V_new)/V_new;
    
    %Do optimization
    %Vopt = CheckOpt(x,R,s_,c,V,PolicyRulesInit,Para);
    %Check2(n) = (Vopt-V_new)/V_new;
    %[Vopt1,Vopt2] = CheckOpt(x,R,s_,c,V,PolicyRules(1:3),Para);
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

%
% % Caption : fig:ValueFunction - This plot depicts the value function
figure()


figure()

for Rctr=1:4
    subplot(2,2,Rctr)
    fplot(@(x) funeval(c(s_,:)',V(s_),[x RList(Rctr)]),[ucbtild_bounds(1) ucbtild_bounds(2)],'-k');
    xlabel('$x$','Interpreter','Latex')
    title(['$\rho=$' num2str(RList(Rctr))],'Interpreter','Latex')
    hold on
end

% % Caption : fig:ValueFunction - This plot depicts the value function
figure()
xlist=linspace(min(Para.xGrid),max(Para.xGrid),4);
for xctr=1:4
    subplot(2,2,xctr)
    fplot(@(R) funeval(c(s_,:)',V(s_),[xlist(xctr) R]),[Rbounds(1) Rbounds(2)],'-k');
    xlabel('$R$','Interpreter','Latex')
    title(['$x=$' num2str(xlist(xctr))],'Interpreter','Latex')
    hold on
end


% % Caption : fig:ValueFunction - This plot depicts the value function
figure()
xlist=0
for xctr=1:1
    fplot(@(R) funeval(c(s_,:)',V(s_),[xlist(xctr) R]),[Rbounds(1) Rbounds(2)],'-k');
    xlabel('$\rho$','Interpreter','Latex')
    title(['$x=$' num2str(xlist(xctr))],'Interpreter','Latex')
    hold on
end
%% Policy Rules entire state space
% Caption : fig:PolicyRules - This plot depicts the $\tilde{b}'_2$ as a function of $\tilde{b}_2$
figxePrime =figure('Name','x');
figEDeltaX=figure('Name','EDeltaX');
figEDeltaR=figure('Name','EDeltaRho');
figVDeltaX=figure('Name','VDeltaX');
figVDeltaR=figure('Name','VDeltaRho');

figBtildePrime =figure('Name','btild');
figRprime=figure('Name','$\rho$');
figFOCRes =figure('Name','FOCRes');
figRR = figure('Name', 'RRGraph');
%figHLLH = figure('Name', 'HLLHPlot');
xFineGrid=linspace(ucbtild_bounds(1),ucbtild_bounds(2),35);
RFineGrid=linspace(Rbounds(1),Rbounds(2),35);
RList=linspace(Rbounds(1),Rbounds(2),4);
xList=linspace(ucbtild_bounds(1),ucbtild_bounds(2),4);
s_=1;
for Rctr=1:4
    for xctr=1:length(xFineGrid)
        R=RList(Rctr);
        x=xFineGrid(xctr);
        [PolicyRulesInit]=GetInitialApproxPolicy([x R s_] ,domain,PolicyRulesStore);
        [PolicyRules, V_new,exitflag,fvec]=CheckGradNAG(x,R,s_,c,V,PolicyRulesInit,Para);
        if exitflag==1
            IndxPrint(xctr)=1;
        else
            IndxPrint(xctr)=0;
        end
        
        FOCRes(xctr)=max(abs(fvec));
        xePrime(xctr,:)=PolicyRules(end-S+1:end);
        Rprime(xctr,:)=PolicyRules(end-2*S+1:end-S);
        BtildePrime(xctr,:)=PolicyRules(end-3*S+1:end-2*S);
        EDeltaX(xctr)=sum(Para.P(s_,:).*xePrime(xctr,:))-x;
        VDeltaX(xctr)=sum(Para.P(s_,:).*(xePrime(xctr,:)-x*ones(1,S)).^2)-EDeltaX(xctr).^2;
        EDeltaR(xctr)=sum(Para.P(s_,:).*Rprime(xctr,:))-R;
        VDeltaR(xctr)=sum(Para.P(s_,:).*(Rprime(xctr,:)-R*ones(1,S)).^2)-EDeltaR(xctr).^2;
    end
    
    
    figure(figEDeltaX)
    subplot(2,2,Rctr)
    plot(xFineGrid, EDeltaX,'k','LineWidth',2)
    xlabel('$x$','Interpreter','Latex')
    ylabel('EDeltaX','Interpreter','Latex')
    title(['$\rho=$' num2str(RList(Rctr))],'Interpreter','Latex')
    
    
    
    
    
    figure(figVDeltaX)
    subplot(2,2,Rctr)
    plot(xFineGrid, VDeltaX,'k','LineWidth',2)
    xlabel('$x$','Interpreter','Latex')
    ylabel('VDeltaX','Interpreter','Latex')
    title(['$\rho=$' num2str(RList(Rctr))],'Interpreter','Latex')
    
    
    
    figure(figEDeltaR)
    subplot(2,2,Rctr)
    plot(xFineGrid, EDeltaR,'k','LineWidth',2)
    xlabel('$x$','Interpreter','Latex')
    ylabel('EDeltaR','Interpreter','Latex')
    title(['$\rho=$' num2str(RList(Rctr))],'Interpreter','Latex')
    
    
    
    figure(figVDeltaR)
    subplot(2,2,Rctr)
    plot(xFineGrid, VDeltaR,'k','LineWidth',2)
    xlabel('$x$','Interpreter','Latex')
    ylabel('VDeltaR','Interpreter','Latex')
    title(['$\rho=$' num2str(RList(Rctr))],'Interpreter','Latex')
    
    
    
    figure(figFOCRes)
    subplot(2,2,Rctr)
    plot(xFineGrid, FOCRes,'k','LineWidth',2)
    if Rctr==1
        legend('g_l','g_h')
    end
    xlabel('$x$','Interpreter','Latex')
    ylabel('FOCRes','Interpreter','Latex')
    title(['$\rho=$' num2str(RList(Rctr))],'Interpreter','Latex')
    
    figure(figBtildePrime)
    subplot(2,2,Rctr)
    for s=1:S
    plot(xFineGrid(logical(IndxPrint)), BtildePrime(logical(IndxPrint),s),'LineWidth',2)
     hold on
    end    
    if Rctr==1
        legend(gShockNames)
    end
    
    hold on
    xlabel('$x$','Interpreter','Latex')
    ylabel('$\tilde{b}_2$','Interpreter','Latex')
    title(['$\rho=$' num2str(RList(Rctr))],'Interpreter','Latex')
    
    figure(figxePrime)
    subplot(2,2,Rctr)
    
    plot(xFineGrid(logical(IndxPrint)), xePrime(logical(IndxPrint),:)- repmat(xFineGrid(logical(IndxPrint))',1,S),'LineWidth',2)
     hold on
        if Rctr==1
        legend(gShockNames)
    end
    hold on
    xlabel('$x$','Interpreter','Latex')
    ylabel('$x(s)-x$','Interpreter','Latex')
    title(['$\rho=$' num2str(RList(Rctr))],'Interpreter','Latex')
    
    figure(figRprime)
    subplot(2,2,Rctr)
    
    plot(xFineGrid(logical(IndxPrint)), Rprime(logical(IndxPrint),:)-RList(Rctr),'LineWidth',2);
    hold on
       
    xlabel('$x$','Interpreter','Latex')
    ylabel('$\rho(s)-\rho$','Interpreter','Latex')
    title(['$\rho=$' num2str(RList(Rctr))],'Interpreter','Latex')
end


for xctr=1:4
    for Rctr=1:length(RFineGrid)
        R=RFineGrid(Rctr);
        x=xList(xctr);
        [PolicyRulesInit]=GetInitialApproxPolicy([x R s_] ,domain,PolicyRulesStore);
        [PolicyRules, V_new,exitflag,fvec]=CheckGradNAG(x,R,s_,c,V,PolicyRulesInit,Para);
        if exitflag==1
            IndxPrint(xctr)=1;
        else
            IndxPrint(xctr)=0;
        end
        
         FOCRes(Rctr)=max(abs(fvec));
        xePrime(Rctr,:)=PolicyRules(end-S+1:end);
        Rprime(Rctr,:)=PolicyRules(end-2*S+1:end-S);
        BtildePrime(Rctr,:)=PolicyRules(end-3*S+1:end-2*S);
        EDeltaX(Rctr)=sum(Para.P(s_,:).*xePrime(Rctr,:))-x;
        VDeltaX(Rctr)=sum(Para.P(s_,:).*(xePrime(Rctr,:)-x*ones(1,S)).^2)-EDeltaX(Rctr).^2;
        EDeltaR(Rctr)=sum(Para.P(s_,:).*Rprime(Rctr,:));
        VDeltaR(Rctr)=sum(Para.P(s_,:).*(Rprime(Rctr,:)-R*ones(1,S)).^2)-EDeltaR(Rctr).^2;
        
    end
    figure(figRR)
    subplot(2,2,xctr);
    
    plot(RFineGrid,Rprime(:,:)-repmat(RFineGrid',1,S),'LineWidth',2);
    hold on
        xlabel('$\rho$','Interpreter','Latex')
    ylabel('$\rho(s)-\rho$','Interpreter','Latex')
    title(['$x=$' num2str(xList(xctr))],'Interpreter','Latex')
end
% 
% print(figFOCRes,'-dpng',[plotpath 'FOCResFullDomain.png'])
% print(figEDeltaX,'-dpng',[plotpath 'EDeltaX.png'])
% print(figEDeltaR,'-dpng',[plotpath 'EDeltaR.png'])
% print(figVDeltaX,'-dpng',[plotpath 'VDeltaX.png'])
% print(figVDeltaR,'-dpng',[plotpath 'VDeltaR.png'])
% print(figxePrime,'-dpng',[plotpath 'xePrimeFullDomain.png'])
% print(figBtildePrime,'-dpng',[plotpath 'BtildePrimeFullDomain.png'])
% print(figRR,'-dpng',[plotpath 'RPrimeFullDomain.png'])
