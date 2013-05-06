function [xhat,Coeff_xhat,Rhat,Coeff_Rhat]=ApproximatePolicyRules(Para,domain,c,PolicyRulesStore,V,trunc)
% the arguments are stored in the coeff.mat files and the last argument
% truc controls the truncation of grid at the end points

% -------------------------------------------------------------------------
% This script uses the policy rules to create non linear approximation of LOM
% of state variables  - x,R
%--------------------------------------------------------------------------


S=Para.sSize;
if nargin==5
    trunc=.9;
end

% Setup the approximation domain
xhatMin=Para.xMin*trunc;
xhatMax=Para.xMax*trunc;
RhatMin=Para.RMin*(2-trunc);
RhatMax=Para.RMax*trunc;
xhatGridSize=Para.xGridSize*2;
RhatGridSize=Para.xGridSize*2;
for s=1:S
xhat(s) = fundefn(Para.ApproxMethod,[Para.OrderOfAppx_x Para.OrderOfApprx_R ] ,[xhatMin RhatMin],[xhatMax RhatMax]);
Rhat(s) = fundefn(Para.ApproxMethod,[Para.OrderOfAppx_x Para.OrderOfApprx_R ] ,[xhatMin RhatMin],[xhatMax RhatMax]);
end

xhatGrid=linspace(xhatMin,xhatMax,xhatGridSize);
RhatGrid=linspace(RhatMin,RhatMax,RhatGridSize);
s_=1;

% Compute the fitted poits
ctr=1;
for xind=1:xhatGridSize
    for Rind=1:RhatGridSize
                R=RhatGrid(Rind);
        x=xhatGrid(xind);
        ApproximationDomain(ctr,:)=[xhatGrid(xind) RhatGrid(Rind)];
        [PolicyRulesInit]=GetInitialApproxPolicy([x R s_] ,domain,PolicyRulesStore);
        [PolicyRules, V_new,exitflag,fvec]=CheckGradNAG(x,R,s_,c,V,PolicyRulesInit,Para);
        if exitflag==1
            IndxPrint(ctr)=1;
        else
            IndxPrint(ctr)=0;
        end
         Rhatprime(ctr,:)=PolicyRules(end-3:end-2);
    xhatprime(ctr,:)=PolicyRules(end-1:end);
    
        
ctr=ctr+1;
    end
end


for s=1:S
        Coeff_xhat(s,:)=funfitxy(xhat(s),ApproximationDomain, xhatprime(:,s));
        Coeff_Rhat(s,:)=funfitxy(Rhat(s),ApproximationDomain, Rhatprime(:,s));
end
end