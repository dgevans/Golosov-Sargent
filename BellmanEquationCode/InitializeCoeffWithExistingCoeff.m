function [ domain, c, PolicyRulesStore] = InitializeCoeffWithExistingCoeff( Para, V,BellmanData)
%% INITIALIZE THE COEFF
% This section uses the stationary policies to initialze the value
% functions. The block prdoduces three outcomes 
% 1. domain which is the vectorized domain
% 2. PolicyRulesStore : This serves as a matrix of guess for optimal
% policies at the interpolation nodes
%3. c0 : initial coeffecients
S = length(Para.P);
xGrid=Para.xGrid;
RGrid=Para.RGrid;
% initialize using deterministic case
Para.g=mean(Para.g)*ones(1,length(Para.g));
Para.theta_1=mean(Para.theta_1)*ones(1,length(Para.theta_1));
Para.theta_2=mean(Para.theta_2)*ones(1,length(Para.theta_2));
lastWokedGuess=0.5*ones(1,2*S-1);
PolicyRulesStore = [];
for s_=1:S
    n=1;
    if true               
        for xctr=1:Para.xGridSize
            for Rctr=1:Para.RGridSize   
                x_=xGrid(xctr);
                R_=RGrid(Rctr);
                %if R_>Rbar(u2btildctr)
                domain_(s_,n,:)=[x_ R_ ];
                % Solve for  c1
                 V0(s_,n) = funeval(BellmanData.c(s_,:)',BellmanData.V(s_),[x_ R_]);
                 [PolicyRulesInit]=GetInitialApproxPolicy([x_ R_ s_] ,BellmanData.domain,BellmanData.PolicyRulesStore);
                 [PolicyRules, V_new,exitflag,fvec]=CheckGradNAG(x_,R_,s_,BellmanData.c,BellmanData.V,PolicyRulesInit,Para);
           
                PolicyRulesStore=vertcat(PolicyRulesStore,PolicyRules) ;
                n=n+1;
%end
            end
        end
        c0(s_,:)=funfitxy(V(s_),squeeze(domain_(s_,:,:)),squeeze(V0(s_,:))' );
    else
        c0(s_,:)=c0(1,:);
        V0(s_,:)=V0(1,:);
    end
end
domain = [];

for s_ = 1:S
    domain =[domain;
             squeeze(domain_(1,:,:)) s_*ones(length(domain_),1)];
end
c=c0;
save([ Para.datapath 'c1.mat' ] , 'c');

end

