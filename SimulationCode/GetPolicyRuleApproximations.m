function [PolicyFunctions,errorPolicyFunctions]=GetPolicyRuleApproximations(Para,c,V,xGridSize,RGridSize,PolicyRulesStoreOld,ApproxMethod,OrderOfAppx_x,OrderOfApprx_R,order,domainOld)
n1=Para.n1;
n2=Para.n2;
alpha_1=Para.alpha_1;
alpha_2=Para.alpha_2;
theta_1=Para.theta_1;
theta_2=Para.theta_2;
psi=Para.psi;
beta=Para.beta;
S=length(Para.P);
xMin=Para.xMin;
xMax=Para.xMax;
RMin=Para.RMin;
RMax=Para.RMax;
domain=cartprod(linspace(xMin,xMax,xGridSize),linspace(RMin,RMax,RGridSize),(1:S));

% Define the functional space
for s_=1:S
for s=1:S
xhat(s_,s)=fundefn(ApproxMethod,[OrderOfAppx_x OrderOfApprx_R ] ,[xMin RMin],[xMax RMax],order);
Rhat(s_,s)=fundefn(ApproxMethod,[OrderOfAppx_x OrderOfApprx_R ] ,[xMin RMin],[xMax RMax],order);
c1hat(s_,s)=fundefn(ApproxMethod,[OrderOfAppx_x OrderOfApprx_R ] ,[xMin RMin],[xMax RMax],order);
c2hat(s_,s)=fundefn(ApproxMethod,[OrderOfAppx_x OrderOfApprx_R ] ,[xMin RMin],[xMax RMax],order);
l1hat(s_,s)=fundefn(ApproxMethod,[OrderOfAppx_x OrderOfApprx_R ] ,[xMin RMin],[xMax RMax],order);
l2hat(s_,s)=fundefn(ApproxMethod,[OrderOfAppx_x OrderOfApprx_R ] ,[xMin RMin],[xMax RMax],order);
end
end

err=[];
try
    matlabpool('size')
catch err
end

if isempty(err)
   
    if(matlabpool('size') == 0)
        matlabpool open local;
        
    end
    
    
end


x_slice=domain(:,1) ;
R_slice=domain(:,2) ;
s_slice=domain(:,3) ;
GridSize=length(domain);
S = length(Para.P);
 
   parfor ctr=1:GridSize
        x=x_slice(ctr) ;
        R=R_slice(ctr) ;
        s_=s_slice(ctr);
        [PolicyRulesInit]=GetInitialApproxPolicy([x R s_] ,domainOld,PolicyRulesStoreOld);

        % INNER OPTIMIZATION
        %[PolicyRules, V_new,exitflag,~]=CheckGradNAG2Shocks(x,R,s_,c,V,xInit',Para);
        [PolicyRulesStore(ctr,:), ~,ExitFlag(ctr,:),~]=CheckGradNAG(x,R,s_,c,V,PolicyRulesInit',Para);
         
    
    end

    for s_=1:S
        IndxUnSolved_s=find(~(ExitFlag==1));
        IndxSolved_s=(s_-1)*GridSize/S+find(ExitFlag((s_-1)*GridSize/S+1:s_*GridSize/S)==1);
        for s=1:S
        Coeff_c1(s_,s,:)=funfitxy(c1hat(s_,s),domain(IndxSolved_s,1:2), PolicyRulesStore(IndxSolved_s,s));
        c1{s_,s}=@(x,R) funeval(squeeze(Coeff_c1(s_,s,:)),c1hat(s),[x R]);
        error_c1{s_,s}=funeval(squeeze(Coeff_c1(s_,s,:)),c1hat(s),domain(IndxSolved_s,1:2))-PolicyRulesStore(IndxSolved_s,s);
        
        Coeff_c2(s_,s,:)=funfitxy(c2hat(s_,s),domain(IndxSolved_s,1:2), PolicyRulesStore(IndxSolved_s,S+s));
        c2{s_,s}=@(x,R) funeval(squeeze(Coeff_c2(s_,s,:)),c2hat(s),[x R]);
        error_c2{s_,s}=funeval(squeeze(Coeff_c2(s_,s,:)),c2hat(s),domain(IndxSolved_s,1:2))-PolicyRulesStore(IndxSolved_s,S+s);
        
        Coeff_l1(s_,s,:)=funfitxy(l1hat(s_,s),domain(IndxSolved_s,1:2), PolicyRulesStore(IndxSolved_s,2*S+s));
        l1{s_,s}=@(x,R) funeval(squeeze(Coeff_l1(s_,s,:)),l1hat(s),[x R]);
        error_l1{s_,s}=funeval(squeeze(Coeff_l1(s_,s,:)),l1hat(s),domain(IndxSolved_s,1:2))-PolicyRulesStore(IndxSolved_s,2*S+s);
        
        Coeff_l2(s_,s,:)=funfitxy(l2hat(s_,s),domain(IndxSolved_s,1:2), PolicyRulesStore(IndxSolved_s,3*S+s));
        l2{s_,s}=@(x,R) funeval(squeeze(Coeff_l2(s_,s,:)),l2hat(s),[x R]);
        error_l2{s_,s}=funeval(squeeze(Coeff_l2(s_,s,:)),l2hat(s),domain(IndxSolved_s,1:2))-PolicyRulesStore(IndxSolved_s,3*S+s);
        
        Coeff_xhat(s_,s,:) = funfitxy(xhat(s_,s),domain(IndxSolved_s,1:2), PolicyRulesStore(IndxSolved_s,end-S+s));
        xPrime{s_,s}=@(x,R) funeval(squeeze(Coeff_xhat(s_,s,:)),xhat(s),[x R]);
        error_xPrime{s_,s}=funeval(squeeze(Coeff_xhat(s_,s,:)),xhat(s),domain(IndxSolved_s,1:2))-PolicyRulesStore(IndxSolved_s,end-S+s);
        
        Coeff_Rhat(s_,s,:) = funfitxy(Rhat(s_,s),domain(IndxSolved_s,1:2),PolicyRulesStore(IndxSolved_s,end-2*S+s));    
        RPrime{s_,s}=@(x,R) funeval(squeeze(Coeff_Rhat(s_,s,:)),Rhat(s),[x R]);
        error_RPrime{s_,s}=funeval(squeeze(Coeff_Rhat(s_,s,:)),Rhat(s),domain(IndxSolved_s,1:2))-PolicyRulesStore(IndxSolved_s,end-2*S+s);
        
        
        end
   end




PolicyFunctions.c1=c1;
PolicyFunctions.c2=c2;
PolicyFunctions.l1=l1;
PolicyFunctions.l2=l2;
PolicyFunctions.xPrime=xPrime;
PolicyFunctions.RPrime=RPrime;


errorPolicyFunctions.c1=error_c1;
errorPolicyFunctions.c2=error_c2;
errorPolicyFunctions.l1=error_l1;
errorPolicyFunctions.l2=error_l2;
errorPolicyFunctions.xPrime=error_xPrime;
errorPolicyFunctions.RPrime=error_RPrime;


end