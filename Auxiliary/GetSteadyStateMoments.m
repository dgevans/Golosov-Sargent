function [res]=GetSteadyStateMoments(Para,c10guess,c20guess)
%% Define the funtional space
V(1) = fundefn(Para.ApproxMethod,[Para.OrderOfAppx_u2btild Para.OrderOfApprx_R ] ,[Para.u2btildMin Para.RMin],[Para.u2btildMax Para.RMax]);
V(2) = V(1);
    %% INITIALIZE THE COEFF
    %  This function computes c1,c2,l1,l2 and the value for an arbitrary x, R.
    % This section solves for V i.e the value function at the end of period 1
    % with g_t=g for all t >1. since the value function is static we need to
    % solve a equation in c_1 for each x,R. Th function getValueC1 does the job
    % by solving for the two roots of this equation and using the one that
    % supports the highest utility
    tic
    gTrue=Para.g;
    Para.g=mean(gTrue)*ones(2,1);
    for s_=1:Para.sSize
        n=1;
        if s_==1
            
            
            for u2btildctr=1:Para.u2btildGridSize
                for Rctr=1:Para.RGridSize
                    
                    u2btild_=Para.u2btildGrid(u2btildctr);
                    R_=Para.RGrid(Rctr);
                    %if R_>Rbar(u2btildctr)
                    x_state_(s_,n,:)=[u2btild_ R_ ];
                    % Solve for  c1
                    
                    c1_=max(getValueC1(u2btild_,R_,s_,Para ),.0001);
                    
                  
                    if c1_<.001
                        ExitFlagT(n)=0;
                    else
                        ExitFlagT(n)=1;
                    end
                    % compute c2
                    c2_=R_^(-1)*c1_;
                    TotalResources=(c1_*Para.n1+c2_*Para.n2+Para.g(s_));
                    FF=R_*Para.theta_2/Para.theta_1;
                    DenL2=Para.n1*Para.theta_1*FF+Para.theta_2*Para.n2;
                    l2_=(TotalResources-Para.n1*Para.theta_1+Para.n1*Para.theta_1*FF)/(DenL2);
                    l1_= 1-FF*(1-l2_);
                    u2btildPrime_=u2btild_;
                    V0(s_,n)=(Para.alpha_1*uBGP(c1_,l1_,Para.psi)+Para.alpha_2*uBGP(c2_,l2_,Para.psi))/(1-Para.beta);
                
                    xInit_0(s_,n,:)=[c1_ c2_ l1_ l2_ u2btildPrime_/(Para.psi*c2_^(-1)) R_ u2btildPrime_];
                    n=n+1;
                    %end
                end
            end
            c0(s_,:)=funfitxy(V(s_),squeeze(x_state_(s_,logical(ExitFlagT==1),:)),squeeze(V0(s_,logical(ExitFlagT==1)))' );
        else
            c0(s_,:)=c0(s_-1,:);
            V0(s_,:)=V0(s_-1,:);
            xInit_0(s_,:)=xInit_0(s_-1,:);
        end
    end
    disp('Number of points solved in initialization')
    sum(ExitFlagT)
    disp('Number of points solved out of a total of ')
    length(ExitFlagT)
    
    Para.g=gTrue;
    x_state=vertcat([squeeze(x_state_(1,:,:)) ones(length(x_state_),1)] ,[squeeze(x_state_(1,:,:)) 2*ones(length(x_state_),1)]);
    scatter(x_state(:,1),x_state(:,2))
    c=c0;
 
    % slicing the state space for parfor loop later
    u2btild_slice=x_state(:,1) ;
    R_slice=x_state(:,2) ;
    s_slice=x_state(:,3) ;
    % This stores the values of the policy functions and multipliers that last
    % worked
    
    PolicyRulesWorked=[xInit_0(1,1,1) xInit_0(2,1,1) xInit_0(1,1,2)];
    
    % This stores the policy rules for each point in the state
    % space.
    PolicyRulesStore1=[squeeze(xInit_0(1,:,1))' squeeze(xInit_0(1,:,1))' ...
        squeeze(xInit_0(1,:,2))' squeeze(xInit_0(1,:,2))'...
        squeeze(xInit_0(1,:,3))' squeeze(xInit_0(1,:,3))' ...
        squeeze(xInit_0(1,:,4))' squeeze(xInit_0(1,:,4))' ....
        squeeze(xInit_0(1,:,5))' squeeze(xInit_0(1,:,5))' ....
        squeeze(xInit_0(1,:,6))' squeeze(xInit_0(1,:,6))' ....
        squeeze(xInit_0(1,:,7))' squeeze(xInit_0(1,:,7))' ....
        ];
    PolicyRulesStore2=[squeeze(xInit_0(2,:,1))' squeeze(xInit_0(2,:,1))' ...
        squeeze(xInit_0(2,:,2))' squeeze(xInit_0(2,:,2))'...
        squeeze(xInit_0(2,:,3))' squeeze(xInit_0(2,:,3))' ...
        squeeze(xInit_0(2,:,4))' squeeze(xInit_0(2,:,4))' ....
        squeeze(xInit_0(2,:,5))' squeeze(xInit_0(2,:,5))' ....
        squeeze(xInit_0(2,:,6))' squeeze(xInit_0(2,:,6))' ....
        squeeze(xInit_0(2,:,7))' squeeze(xInit_0(2,:,7))' ....
        ];
    PolicyRulesStore=vertcat(PolicyRulesStore1,PolicyRulesStore2);












% SOLVE THE T-0 PROBLEM given btild(-1)
btild_1=Para.btild_1;
disp('Computed V...Now solving V0(btild_1) where btild_1 is')
disp(btild_1)
% c1 and c2 solve
options=optimset('Display','off');
[x,~,exitflagv0,~,~] = fminunc(@(x)  getValue0(x, btild_1,1,Para,c,V),[ c10guess c20guess],options);
if ~(exitflagv0==1)
[x,~,exitflagv0,~,~] = fminunc(@(x)  getValue0(x, btild_1,1,Para,c,V),[ 1 1/Para.RMax],options);
end

if ~(exitflagv0==1)
    disp('Optimization failed for V0 once ..trying with fmincon')
    opts = optimset('Algorithm', 'interior-point', 'Display','off', ...
        'GradObj','off','GradConstr','off',...
        'MaxIter',1000, ...
        'TolX', Para.ctol/10, 'TolFun', Para.ctol, 'TolCon', Para.ctol,'MaxTime',200);
    lb=[0.001 0.001];
    ub=[10 10];
    %[x,fval,exitflagv0,output,lambda]  =fmincon(@(x) getValue0(x, btild_1,1,Para,c,V),[ x ],[],[],[],[],lb,ub,[],opts);
    [x,~,exitflagv0,output,lambda]  =ktrlink(@(x) getValue0(x, btild_1,1,Para,c,V),[ c10guess c20guess],[],[],[],[],lb,ub,[],opts);
    
end
c10 = x(1);
c20 = x(2);
R0=c10/c20;
TotalResources=(c10*Para.n1+c20*Para.n2+Para.g(1));
FF=R0*Para.theta_2/Para.theta_1;
DenL2=Para.n1*Para.theta_1*FF+Para.theta_2*Para.n2;
l20=(TotalResources-Para.n1*Para.theta_1+Para.n1*Para.theta_1*FF)/(DenL2);
l10= 1-FF*(1-l20);
BracketTerm=l20/(1-l20)-(l10/(1-l10))*R0;
u2btildprime0=(((1-Para.psi)/(Para.psi))*BracketTerm+btild_1/(Para.beta*Para.psi)+R0-1)*Para.psi;
btildprime0=u2btildprime0/(c20^-1*Para.psi) ;
Rprime0=c20^(-1)/c10^(-1);
ul20=(1-Para.psi)/(1-l20);
ul10=(1-Para.psi)/(1-l10);
uc20=Para.psi/c20;
uc10=Para.psi/c10;
tau0=1-(ul20/(Para.theta_2*uc20));

 u2btild_=u2btildprime0;
                    R_=Rprime0;
                    c1=max(getValueC1(u2btild_,R_,s_,Para ),.0001);
                   
                    % compute c2
                    c2=R_^(-1)*c1;
                    TotalResources=(c1*Para.n1+c2_*Para.n2+Para.g(s_));
                    FF=R_*Para.theta_2/Para.theta_1;
                    DenL2=Para.n1*Para.theta_1*FF+Para.theta_2*Para.n2;
                    l2=(TotalResources-Para.n1*Para.theta_1+Para.n1*Para.theta_1*FF)/(DenL2);
                    l1= 1-FF*(1-l2);
                   ul2=(1-Para.psi)/(1-l2);
ul1=(1-Para.psi)/(1-l1);
uc1=Para.psi/c1;
uc2=Para.psi/c2;
tau1=1-(ul2/(Para.theta_2*uc2));
 Trans=c2-l2.*ul2./uc2;
    
  
    
     % Income
    %AfterTaxWageIncome_Agent2=l2.*ul2./uc2;
    %AfterTaxWageIncome_Agent1=l1.*ul1./uc1;
    AfterTaxWageIncome_Agent2=l2.*Para.theta_2;
    AfterTaxWageIncome_Agent1=l1.*Para.theta_1;
    % Gini Coeff
    GiniCoeff=(AfterTaxWageIncome_Agent2 +2*AfterTaxWageIncome_Agent1)./(AfterTaxWageIncome_Agent2+AfterTaxWageIncome_Agent1)-3/2;
  
AvgHrs=(l1+l2)/2;
res.tau1=tau1;
res.l1=l1;
res.l2=l2;
res.GiniCoeff=GiniCoeff;
res.AvgHrs=AvgHrs;
res.AvgFrishElasticity=(Para.n1*(1/l1-1)+Para.n2*(1/l2-1))/(Para.n1+Para.n2);
res.exitflagv0=exitflagv0;
res.VarHrs=var([log(l1) log(l2)],1);
res.VarLogWages=var([log(AfterTaxWageIncome_Agent2) log(AfterTaxWageIncome_Agent1)],1);
res.c10=c10;
res.c20=c20;
end