% CAPTION : fig:flagPoints - This plots the sucess of the optimizer to
% solve the FOC at the points selected in the state space for the final set of coeffecients. The red points
% denote failure.
function GetPlotsForFinalSolution(Para)
load([Para.datapath Para.StoreFileName])
close all;
olddatapath=Para.datapath;
oldtexpath=Para.texpath;
oldplotpath=Para.plotpath;
plotpath='Graphs/Calibration/temp/';
datapath='Data/Calibration/';
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


u2btildLL=Para.u2btildMin;
u2btildUL=Para.u2btildMax;
ucbtild_bounds = [u2btildLL,u2btildUL];
Rbounds=[min(Para.RGrid),max(Para.RGrid)];

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
    fplot(@(u2btild) funeval(c(s_,:)',V(s_),[u2btild RList(Rctr)]),[ucbtild_bounds(1) ucbtild_bounds(2)],'-k');
end
xlabel('$x$','Interpreter','Latex')
title(['$R=$' num2str(RList(Rctr))],'Interpreter','Latex')
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
    title(['$R=$' num2str(RList(Rctr))],'Interpreter','Latex')
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
    xlabel('$R$','Interpreter','Latex')
    title(['$x=$' num2str(xlist(xctr))],'Interpreter','Latex')
    hold on
end

print(gcf,'-dpng',[plotpath 'ValueFunctionRx_0.png'])


%% Policy Rules entire state space
% Caption : fig:PolicyRules - This plot depicts the $\tilde{b}'_2$ as a function of $\tilde{b}_2$
figu2BtildePrime =figure('Name','x');
figBtildePrime =figure('Name','btild');
figRprime=figure('Name','R');
figFOCRes =figure('Name','FOCRes');
u2bdiffFineGrid=linspace(min(Para.u2bdiffGrid),max(Para.u2bdiffGrid),35);
RList=linspace(min(Para.RGrid),max(Para.RGrid),4);
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
        Rprime(u2btildctr,:)=PolicyRules(end-3:end-2);
    end
    
    figure(figFOCRes)
    subplot(2,2,Rctr)
    plot(u2bdiffFineGrid, FOCRes,'k')
    if Rctr==1
        legend('g_l','g_h')
    end
    xlabel('$x$','Interpreter','Latex')
    ylabel('FOCRes','Interpreter','Latex')
    title(['$R=$' num2str(RList(Rctr))])
    
    figure(figBtildePrime)
    subplot(2,2,Rctr)
    plot(u2bdiffFineGrid(logical(IndxPrint)), BtildePrime(logical(IndxPrint),1),'k')
    hold on
    plot(u2bdiffFineGrid(logical(IndxPrint)), BtildePrime(logical(IndxPrint),2),':k')
    hold on
    if Rctr==1
        legend('g_l','g_h')
    end
    %plot(u2bdiffFineGrid, 0*u2bdiffFineGrid,':k');
    hold on
    %plot(u2bdiffFineGrid,repmat([u2btildLL u2btildUL],length(u2bdiffFineGrid),1)-[u2bdiffFineGrid' u2bdiffFineGrid'] ,':r')
    %
    xlabel('$x$','Interpreter','Latex')
    ylabel('$\tilde{b}_2$','Interpreter','Latex')
    title(['$R=$' num2str(RList(Rctr))])
    
    figure(figu2BtildePrime)
    subplot(2,2,Rctr)
    plot(u2bdiffFineGrid(logical(IndxPrint)), u2BtildePrime(logical(IndxPrint),1)- u2bdiffFineGrid(logical(IndxPrint))','k')
    hold on
    plot(u2bdiffFineGrid(logical(IndxPrint)), u2BtildePrime(logical(IndxPrint),2)-u2bdiffFineGrid(logical(IndxPrint))',':k')
    hold on
    if Rctr==1
        legend('g_l','g_h')
    end
    %plot(u2bdiffFineGrid, 0*u2bdiffFineGrid,':k');
    hold on
    %plot(u2bdiffFineGrid,repmat([u2btildLL u2btildUL],length(u2bdiffFineGrid),1)-[u2bdiffFineGrid' u2bdiffFineGrid'] ,':r')
    %
    xlabel('$x$','Interpreter','Latex')
    ylabel('$x*-x$','Interpreter','Latex')
    title(['$R=$' num2str(RList(Rctr))])
    
    
    %
    figure(figRprime)
    subplot(2,2,Rctr)
    plot(u2bdiffFineGrid(logical(IndxPrint)), Rprime(logical(IndxPrint),1),'k');
    hold on
    plot(u2bdiffFineGrid(logical(IndxPrint)), Rprime(logical(IndxPrint),2),':k');
    xlabel('$x$','Interpreter','Latex')
    ylabel('$R^{*}$','Interpreter','Latex')
    title(['$R=$' num2str(RList(Rctr))])
end
print(figFOCRes,'-dpng',[plotpath 'FOCResFullDomain.png'])
print(figu2BtildePrime,'-dpng',[plotpath 'u2BtildePrimeFullDomain.png'])
print(figBtildePrime,'-dpng',[plotpath 'BtildePrimeFullDomain.png'])
print(figRprime,'-dpng',[plotpath 'RPrimeFullDomain.png'])




end
