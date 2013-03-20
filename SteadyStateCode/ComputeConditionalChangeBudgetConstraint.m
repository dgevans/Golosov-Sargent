
function NormDecomp=ComputeConditionalChangeBudgetConstraint(Para,c,V,s_,domain,PolicyRulesStore,x,R)
P=Para.P;
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
xLL=Para.xMin;
xUL=Para.xMax;

        [PolicyRulesInit]=GetInitialApproxPolicy([x R s_] ,domain,PolicyRulesStore);
        [PolicyRules, V_new,exitflag,fvec]=CheckGradNAG(x,R,s_,c,V,PolicyRulesInit,Para);            
    %PolicyRules=[c1_1 c1_2 c2_1 c2_2 l1(1) l1(end) l2(1) l2(end) debtprime c2_1^(-1)/c1_1^(-1) c2_2^(-1)/c1_2^(-1) xprime(1) xprime(end)]
    c1=PolicyRules(1:S);
    c2=PolicyRules(S+1:2*S);
    l1=PolicyRules(2*S+1:3*S);
    l2=PolicyRules(3*S+1:4*S);
    ul2=(1-psi)./(1-l2);
    uc2=psi./(c2.^(sigma));
    ul1=(1-psi)./(1-l1);
    uc1=psi./(c1.^(sigma));
    Rprime=PolicyRules(end-2*S+1:end-S);
    % x' - u_c_2* debtprime
    xprime=PolicyRules(end-S+1:end);
    % debt 
    debtprime=-PolicyRules(end-3*S+1:end-2*S);
    
    DeltaLogInt=abs(log(uc2(1))-log(uc2(S)))*100;
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
    DebtDiffHist=n1*(debtprime(end)-debtprime(1));

        ETau=Para.P(s_,:)*Tau';
        SigmaTau=(Para.P(s_,:)*Tau.^2'-ETau^2)^.5;
        
        ETrans=Para.P(s_,:)*(Trans./y)';
        SigmaTrans=(Para.P(s_,:)*(Trans./y).^2'-ETrans^2)^.5;
        
       
        Edebtprime=Para.P(s_,:)*(debtprime./y)';
        Sigmadebtprime=(Para.P(s_,:)*(debtprime./y).^2'-Edebtprime^2)^.5;
        
%     GiniCoeff=(AfterTaxWageIncome_Agent2 +2*AfterTaxWageIncome_Agent1)./(AfterTaxWageIncome_Agent2+AfterTaxWageIncome_Agent1)-3/2;
             LaborTaxRates=Tau;
TransfersGDPRatio=Trans*(n1+n2);
DebtGDPRatio =debtprime.*Para.n1/y;
GDPDecompostion=[n1*sum(Para.P(s_,:).*c1./y) n2*sum(Para.P(s_,:).*c2./y) sum(Para.P(s_,:).*g./y)];

        
        
L2=(alpha_2*g + alpha_1*theta_2 - alpha_2*theta_1 - alpha_2*g*psi + alpha_2*psi*theta_1 + alpha_2*psi*theta_2)/(alpha_1*theta_2 + alpha_2*theta_2);
L1=1-(theta_2/theta_1)*(alpha_1/alpha_2)*(1-L2);
C2=(psi/(1-psi))*theta_2*(1-L2);
C1=(alpha_1/alpha_2)*C2;


Y=theta_1*L1+theta_2*L2;
YBar=sum(P(s_,:).*Y);
NormDecomp=-([Para.g(end)-Para.g(1),n1*(-debtprime(end)+debtprime(1)),(n1+n2)*(Trans(end)-Trans(1)),n1*(theta_1(end).*Tau(end).*l1(end)- theta_1(1).*Tau(1).*l1(1)),Para.n2*(theta_2(end).*Tau(end).*l2(end)- theta_2(1).*Tau(1).*l2(1)) y(end)-y(1)]/YBar)*100;
