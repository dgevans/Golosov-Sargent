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
% initialize using deterministic case
Para.g=mean(Para.g)*ones(1,length(Para.g));
Para.theta_1=mean(Para.theta_1)*ones(1,length(Para.theta_1));
Para.theta_2=mean(Para.theta_2)*ones(1,length(Para.theta_2));
Para.beta=mean(Para.beta)*ones(1,length(Para.beta));
lastWokedGuess=0.5*ones(1,2*S-1);
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
                cRat = R_^(-1/Para.sigma);
                c1 = (0.8*(Para.n1*Para.theta_1+Para.n2*Para.theta_2)-Para.g)./(Para.n1+cRat*Para.n2);
                c1 = c1(:)';
                c2_ = cRat*c1; c2_(S) = []; 
                options = optimset('Display','off');
                warning off
                [xSS,~,exitFlag] = fsolve(@(z) SteadyStateResiduals(z,x_,R_,Para,s_),[c1 c2_],options);
                [res, c1, c2, l1, l2] = SteadyStateResiduals(xSS,x_,R_,Para,s_);
%                 if(exitFlag ~= 1)
%                  [xSS,~,exitFlag] = fsolve(@(z) SteadyStateResiduals(z,x_,R_,Para,s_),[lastWokedGuess],options);
%                 [res, c1, c2, l1, l2] = SteadyStateResiduals(xSS,x_,R_,Para,s_);
%                 end
                if(exitFlag ~= 1)
                    R_
                    x_
                    res
               %     throw(MException('Could Not Find Steady State Policies'));
                else
                    lastWokedGuess=[c1(:)' c2(1:S-1)];
                end
                
                V0(s_,n) = (Para.alpha_1*uAlt(c1,l1,Para.psi,Para.sigma)+Para.alpha_2*uAlt(c2,l2,Para.psi,Para.sigma)./(1-Para.beta))*Para.P(s_,:)';
                
                
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
domain = [];
PolicyRulesStore = [];
for s_ = 1:S
    domain =[domain;
             squeeze(domain_(1,:,:)) s_*ones(length(domain_),1)];
    PolicyRulesStore = [PolicyRulesStore;
              squeeze(xInit_0(s_,:,:))];
end
c=c0;
save([ Para.datapath 'c1.mat' ] , 'c');

end

