% CAPTION : fig:flagPoints - This plots the sucess of the optimizer to
% solve the FOC at the points selected in the state space for the final set of coeffecients. The red points
% denote failure.
function ComputeConditionalMoments(Para,Domain)
flagPlot2PeriodDrifts=0;
oldplotpath=Para.plotpath;
load([Para.datapath Para.StoreFileName])
close all;
plotpath=oldplotpath;
mkdir(plotpath);
datapath='Data/Calibration/';
Para.P
n1=Para.n1;
n2=Para.n2;
alpha_1=Para.alpha_1;
alpha_2=Para.alpha_2;
%disp('Govt Exp')
g=Para.g
theta_1=Para.theta_1;
theta_2=Para.theta_2;
psi=Para.psi;
beta=Para.beta;
sigma=Para.sigma;
S=length(Para.P);
xSolved=domain(IndxSolved,:);
xUnSolved=domain(IndxUnSolved,:);
xLL=Para.xMin;
xUL=Para.xMax;
ucbtild_bounds = [xLL,xUL]*.9;
Rbounds=[min(Para.RGrid),max(Para.RGrid)];

%% Policy Rules entire state space
% Caption : fig:PolicyRules - This plot depicts the $\tilde{b}'_2$ as a
% function of $\tilde{b}_2$
figETau=figure();
figSigmaTau=figure();
figETrans=figure();
figSigmaTrans=figure();
figSigmaBtildePrime=figure();
figEBtildePrime=figure();
figDecomp=figure();
figGDPDecomp=figure();
figDecomp_y=figure();
xFineGrid=linspace(ucbtild_bounds(1),ucbtild_bounds(end),35);
RFineGrid=linspace(Rbounds(1),Rbounds(end),35);
RList=linspace(Rbounds(1),Rbounds(end),4);
xList=linspace(ucbtild_bounds(1),ucbtild_bounds(end),4);
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
        
            
            %PolicyRules=[c1_1 c1_2 c2_1 c2_2 l1(1) l1(end) l2(1) l2(end) btildprime c2_1^(-1)/c1_1^(-1) c2_2^(-1)/c1_2^(-1) xprime(1) xprime(end)]
    c1=PolicyRules(1:S);
    c2=PolicyRules(S+1:2*S);
    l1=PolicyRules(2*S+1:3*S);
    l2=PolicyRules(3*S+1:4*S);
    ul2=(1-psi)./(1-l2);
    uc2=psi./(c2.^(sigma));
    ul1=(1-psi)./(1-l1);
    uc1=psi./(c1.^(sigma));
    Rprime=PolicyRules(end-2*S+1:end-S);
    % x' - u_c_2* btildprime
    xprime=PolicyRules(end-S+1:end);
    % btildprime - x'/u_c2
    btildprime=-PolicyRules(end-3*S+1:end-2*S);
    
    
    % TAU - From the WAGE optimality of Agent 2
    Tau=1-(ul1./(theta_1.*uc1));
    
    % OUTPUT
    
    y=c1*n1+c2*n2+g;
    
    
    % TRANSFERS
    % These are transfers computed on the assumption that Agent 2 cannot
    % borrow and lend. The transfers are the difference between his
    % consumption and after tax earning (l . U_l/U_c)
    Trans=c2-l2.*ul2./uc2;
    
    
     % Income
    AfterTaxWageIncome_Agent2=l2.*ul2./uc2+Trans;
    AfterTaxWageIncome_Agent1=l1.*ul1./uc1+Trans;
    
        TransDiffHist=(Trans(end)-Trans(1));
    % diff in labortax agent1
    LaborTaxAgent1DiffHist=theta_1(end)*l1(1)*Tau(end)*n1 - theta_1(1)*l1(1)*Tau(1)*n1;
    % diff in labortax agent2
    LaborTaxAgent2DiffHist=theta_2(end)*l2(S)*Tau(end)*n2 - theta_2(1)*l2(1)*Tau(1)*n2;
      % diff in borrowing
    DebtDiffHist=n1*(btildprime(end)-btildprime(1));

        ETau(xctr)=Para.P(s_,:)*Tau';
        SigmaTau(xctr)=(Para.P(s_,:)*Tau.^2'-ETau(xctr)^2)^.5;
        
        ETrans(xctr)=Para.P(s_,:)*(Trans./y)';
        SigmaTrans(xctr)=(Para.P(s_,:)*(Trans./y).^2'-ETrans(xctr)^2)^.5;
        
       
        Ebtildprime(xctr)=Para.P(s_,:)*(btildprime./y)';
        Sigmabtildprime(xctr)=(Para.P(s_,:)*(btildprime./y).^2'-Ebtildprime(xctr)^2)^.5;
        
%     GiniCoeff=(AfterTaxWageIncome_Agent2 +2*AfterTaxWageIncome_Agent1)./(AfterTaxWageIncome_Agent2+AfterTaxWageIncome_Agent1)-3/2;
     Decomp(xctr,:)=[Para.g(end)-Para.g(1),n1*(btildprime(end)-btildprime(1)),(n1+n2)*(Trans(end)-Trans(1)),n1*(theta_1(end).*Tau(end).*l1(end)- theta_1(1).*Tau(1).*l1(1)),Para.n2*(theta_2(end).*Tau(end).*l2(end)- theta_2(1).*Tau(1).*l2(1))];
             if length(Para.g)>1 
     DecompRatio(xctr,:)=[n1*(btildprime(end)-btildprime(1)),-(n1+n2)*(Trans(end)-Trans(1)),n1*(theta_1(end).*Tau(end).*l1(end)- theta_1(1).*Tau(1).*l1(1)),Para.n2*(theta_2(end).*Tau(end).*l2(end)- theta_2(1).*Tau(1).*l2(1))]*100./(Para.g(end)-Para.g(1));
             else
                 DecompRatio(xctr,:)=-[-n1*(btildprime(end)-btildprime(1)),-(n1+n2)*(Trans(end)-Trans(1)),n1*(theta_1(end).*Tau(end).*l1(end)- theta_1(1).*Tau(1).*l1(1)),Para.n2*(theta_2(end).*Tau(end).*l2(end)- theta_2(1).*Tau(1).*l2(1))]*100;
 
             end
             DeltaTTD_y(xctr,:)=-[(btildprime(end)/y(end)-btildprime(1)./y(1)),(Trans(end)/y(end)-Trans(1)/y(1)),Tau(end)-Tau(1)]*100;
          
             LaborTaxRates(xctr,:)=Tau;
TransfersGDPRatio(xctr,:)=Trans*(n1+n2);
DebtGDPRatio(xctr,:) =btildprime.*Para.n1/y;
GDPDecompostion(xctr,:)=[n1*sum(Para.P(s_,:).*c1./y) n2*sum(Para.P(s_,:).*c2./y) sum(Para.P(s_,:).*g./y)];

    end
    
    
    figure(figDecomp_y)
    subplot(2,2,Rctr)
plot(xFineGrid, DeltaTTD_y(:,1),'r','LineWidth',2)
hold on
plot(xFineGrid, DeltaTTD_y(:,2),'b','LineWidth',2)
hold on
plot(xFineGrid, DeltaTTD_y(:,3),'k','LineWidth',2)
hold on
plot(xFineGrid, DeltaTTD_y(:,3)*0,':k')
if Rctr==1
legend('\Delta debt-gdp',' \Delta transfers-gdp' ,' \Delta tau')
end
xlabel(['$x$'],'Interpreter','Latex')
title(['$\rho=$' num2str(RList(Rctr))],'Interpreter','Latex')


    
    figure(figDecomp)
    subplot(2,2,Rctr)
area(xFineGrid, DecompRatio)
hold on
if Rctr==1
legend('Debt','Transfers','Agent 1 (Taxes)','Agent 2 (Taxes)')
end
xlabel(['$x$'],'Interpreter','Latex')
title(['$\rho=$' num2str(RList(Rctr))],'Interpreter','Latex')
    

figure(figGDPDecomp)
    subplot(2,2,Rctr)
area(xFineGrid, GDPDecompostion)
hold on
if Rctr==1
legend('Agent 1 Consumption','Agent 2 Consumption','Public Consumption')
end
xlabel(['$x$'],'Interpreter','Latex')
title(['$\rho=$' num2str(RList(Rctr))],'Interpreter','Latex')
    

    

    
    figure(figSigmaTau)
    subplot(2,2,Rctr)
    plot(xFineGrid, SigmaTau,'k','LineWidth',2)
    xlabel('$x$','Interpreter','Latex')
    ylabel('Volatility of Labor Tax','Interpreter','Latex')
    title(['$\rho=$' num2str(RList(Rctr))],'Interpreter','Latex')
    
    figure(figETau)
    subplot(2,2,Rctr)
    plot(xFineGrid, ETau,'k','LineWidth',2)
    xlabel('$x$','Interpreter','Latex')
    ylabel('Average Labor Tax','Interpreter','Latex')
    title(['$\rho=$' num2str(RList(Rctr))],'Interpreter','Latex')
    
    
    
    figure(figSigmaTrans)
    subplot(2,2,Rctr)
    plot(xFineGrid, SigmaTrans,'k','LineWidth',2)
    xlabel('$x$','Interpreter','Latex')
    ylabel('Volatility of transfers-gdp','Interpreter','Latex')
    title(['$\rho=$' num2str(RList(Rctr))],'Interpreter','Latex')
    
    figure(figETrans)
    subplot(2,2,Rctr)
    plot(xFineGrid, ETrans,'k','LineWidth',2)
    xlabel('$x$','Interpreter','Latex')
    ylabel('Average transfers-gdp','Interpreter','Latex')
    title(['$\rho=$' num2str(RList(Rctr))],'Interpreter','Latex')
    
    
    figure(figSigmaBtildePrime)
    subplot(2,2,Rctr)
    plot(xFineGrid, Sigmabtildprime,'k','LineWidth',2)
    xlabel('$x$','Interpreter','Latex')
    ylabel('Volatility of debt-gdp','Interpreter','Latex')
    title(['$\rho=$' num2str(RList(Rctr))],'Interpreter','Latex')
    
    figure(figEBtildePrime)
    subplot(2,2,Rctr)
    plot(xFineGrid, Ebtildprime,'k','LineWidth',2)
    xlabel('$x$','Interpreter','Latex')
    ylabel('Average  debt-gdp','Interpreter','Latex')
    title(['$\rho=$' num2str(RList(Rctr))],'Interpreter','Latex')
    
    
    
end
   print(figETau,'-dpng',[plotpath 'figETau.png']) 
   print(figETrans,'-dpng',[plotpath 'figETrans.png']) 
   print(figSigmaTau,'-dpng',[plotpath 'figSimgaTau.png']) 
   print(figSigmaTrans,'-dpng',[plotpath 'figSimgaTrans.png']) 
   print(figEBtildePrime,'-dpng',[plotpath 'figEBtildePrime.png']) 
   print(figSigmaBtildePrime,'-dpng',[plotpath 'figSigmaBtildePrime.png']) 
   print(figDecomp_y,'-dpng',[plotpath 'figDecomp_y.png']) 
   
 
    