function  GetComparitiveStaticsPlot(Param,ParamName,ParamGrid,ParamBMValue,plotpath)
I=length(Param);
s_=1
for i=1:I
    Para=Param(i)
% load the coeffecients
sigma=Para.sigma;
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

% Compute the SS policies
% Checking against david's code
[ RSS,xSS,PolicyRules ] = findSteadyState( 0,mean(Para.RGrid),Para);

     c1=PolicyRules(1:2);
     c2=PolicyRules(3:4);
     l1=PolicyRules(5:6);
     l2=PolicyRules(7:8);
     ul2=(1-psi)./(1-l2);
     uc2=(psi*c2.^(-sigma));
   Euc2=sum(uc2.*Para.P(s_,:));
     ul1=(1-psi)./(1-l1);
     uc1=(psi*c1.^(-sigma));
     Rprime=[RSS RSS];
%     % x' - u_c_2* btildprime
     u2btildprime=[xSS xSS];
%     % btildprime - x'/u_c2
     btildprime=u2btildprime./uc2;
%
     
%     % TAU - From the WAGE optimality of Agent 2
    Tau=1-(ul2./(theta_2.*uc2));
%     
%     % OUTPUT
     y(1)=c1(1)*n1+c2(1)*n2+g(1);
     y(2)=c1(2)*n1+c2(2)*n2+g(2);
%     
%     % TRANSFERS
%     % These are transfers computed on the assumption that Agent 2 cannot
%     % borrow and lend. The transfers are the difference between his
%     % consumption and after tax earning (l . U_l/U_c)
    Trans=c2-l2.*ul2./uc2;
%     
%      % Income
     AfterTaxWageIncome_Agent2=l2.*ul2./uc2;
     AfterTaxWageIncome_Agent1=l1.*ul1./uc1;
%     % Gini Coeff
Agent2BudgetCheck=c2-Trans-theta_2*l2.*(1-Tau);
Agent1BudgetCheck=c1-Trans-(1-Tau).*theta_1.*l1-btildprime+xSS/(beta*Euc2);
ResourceConsCheck=y-theta_1*n1*l1-theta_2*n2*l2;
GBCCheck=g+Trans*(n1+n2)-xSS/(beta*Euc2)*n1-Tau.*y+btildprime*n1;

%     GiniCoeff=(AfterTaxWageIncome_Agent2 +2*AfterTaxWageIncome_Agent1)./(AfterTaxWageIncome_Agent2+AfterTaxWageIncome_Agent1)-3/2;
     Decomp(i,:)=[Para.g(2)-Para.g(1),n1*(btildprime(2)-btildprime(1)),(n1+n2)*(Trans(2)-Trans(1)),n1*(theta_1.*Tau(2).*l1(2)- theta_1.*Tau(1).*l1(1)),Para.n2*(theta_2.*Tau(2).*l2(2)- theta_2.*Tau(1).*l2(1))];
          DecompRatio(i,:)=[-n1*(btildprime(2)-btildprime(1)),-(n1+n2)*(Trans(2)-Trans(1)),n1*(theta_1.*Tau(2).*l1(2)- theta_1.*Tau(1).*l1(1)),Para.n2*(theta_2.*Tau(2).*l2(2)- theta_2.*Tau(1).*l2(1))]*100./(Para.g(2)-Para.g(1));
LaborTaxRates(i,:)=Tau;
TransfersGDPRatio(i,:)=Trans*(n1+n2);
DebtGDPRatio(i,:) =btildprime.*Para.n1/y;
SteadyState(i,:)=[xSS RSS];
GDPDecompostion(i,:)=[n1*sum(Para.P(s_,:).*c1./y) n2*sum(Para.P(s_,:).*c2./y) sum(Para.P(s_,:).*g./y)];
end
    
figure()
subplot (2,1,1)
area(ParamGrid, DecompRatio)
hold on
  vline(ParamBMValue,'k')

legend('Debt','Transfers','Agent 1 (Taxes)','Agent 2 (Taxes)')
xlabel(['$\' ParamName '$'],'Interpreter','Latex')
ylabel('g(h)-g(l)')

subplot(2,1,2)
area(ParamGrid,GDPDecompostion)
vline(ParamBMValue,'k')
legend('Agent 1 Consumption','Agent 2 Consumption','Public Consumption')
xlabel(['$\' ParamName '$'],'Interpreter','Latex')
ylabel('GDP')
print('-dpng',[ plotpath 'CompStatics' ParamName '.png'])
end

