function [ domain, c, PolicyRulesStore] = InitializeCoeff( Para, V)
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
for s_=1:Para.sSize
    n=1;
    if s_==1                
        for xctr=1:Para.xGridSize
            for Rctr=1:Para.RGridSize   
                x_=xGrid(xctr);
                R_=RGrid(Rctr);
                %if R_>Rbar(u2btildctr)
                domain_(s_,n,:)=[x_ R_ ];
                % Solve for  c1
                cRat = R_^(-1/Para.sigma);
                c1 = (0.8*(Para.n1*Para.theta_1+Para.n2*Para.theta_2)-Para.g)/(Para.n1+cRat*Para.n2);
                c1 = c1(:)';
                c2_ = cRat*c1; c2_(S) = [];
                options = optimset('Display','off');
                [xSS,~,exitFlag] = fsolve(@(z) SteadyStateResiduals(z,x_,R_,Para,s_),[c1 c2_],options);
                [res, c1, c2, l1, l2] = SteadyStateResiduals(xSS,x_,R_,Para,s_);
                if(exitFlag ~= 1)
                    R_
                    x_
                    res
                    throw(MException('Could Not Find Steady State Policies'));
                end
                
                V0(s_,n) = (Para.alpha_1*uAlt(c1,l1,Para.psi,Para.sigma)+Para.alpha_2*uAlt(c2,l2,Para.psi,Para.sigma))*Para.P(s_,:)'/(1-Para.beta);
                
                
                xInit_0(s_,n,:)=[c1 c2 l1 l2 x_*ones(1,S) R_*ones(1,S) x_*ones(1,S)];
                n=n+1;
                %end
            end
        end
        % Initialize the coeffecients via a routine provided
        % by compEcon library  - funfitxy
        c0(s_,:)=funfitxy(V(s_),squeeze(domain_(s_,:,:)),squeeze(V0(s_,:))' );
    else
        c0(s_,:)=c0(1,:);
        V0(s_,:)=V0(1,:);
        xInit_0(s_,:)=xInit_0(1,:);      
    end
end
x_state = [];
PolicyRulesStore = [];
for s_ = 1:S
    domain =[x_state;
             squeeze(domain_(1,:,:)) s_*ones(length(domain_),1)];
    PolicyRulesStore = [PolicyRulesStore;
              squeeze(xInit_0(s_,:,:))];
end
c=c0;
save([ Para.datapath 'c1.mat' ] , 'c');

    
    


end

