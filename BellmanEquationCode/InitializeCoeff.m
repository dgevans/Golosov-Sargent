function [ domain, c, PolicyRulesStore] = InitializeCoeff( Para, V)
%% INITIALIZE THE COEFF
% This section uses the stationary policies to initialze the value
% functions. The block prdoduces three outcomes 
% 1. domain which is the vectorized domain
% 2. PolicyRulesStore : This serves as a matrix of guess for optimal
% policies at the interpolation nodes
%3. c0 : initial coeffecients

xGrid=Para.xGrid;
RGrid=Para.RGrid;
gTrue=Para.g;
Para.g=mean(gTrue)*ones(2,1);
for s_=1:Para.sSize
    n=1;
    if s_==1                
        for xctr=1:Para.xGridSize
            for Rctr=1:Para.RGridSize   
                x_=xGrid(xctr);
                R_=RGrid(Rctr);
                domain_(s_,n,:)=[x_ R_ ];
                % Initialize the guess for stationary policies
                cRat = R_^(-1/Para.sigma);
                c1_1 = (0.8*(Para.n1*Para.theta_1+Para.n2*Para.theta_2)-Para.g(1))/(Para.n1+cRat*Para.n2);
                c1_2 = (0.8*(Para.n1*Para.theta_1+Para.n2*Para.theta_2)-Para.g(2))/(Para.n1+cRat*Para.n2);
                c2_1 = cRat*c1_1;
                options = optimset('Display','off', 'TolFun', 1e-13);
                % Compute the Stationary Plicies using the
                % SteadyStateResiduals routine. 
                [xSS,~,exitFlag] = fsolve(@(x) SteadyStateResiduals(x,x_,R_,Para,s_),[c1_1 c1_2 c2_1],options);
                [res, c1_, c2_, l1_, l2_] = SteadyStateResiduals(xSS,x_,R_,Para,s_);
                % Dispaly points that encoutered errors
                if(exitFlag ~= 1)
                    R_
                    x_
                    res
                end
                % Presented discounted value for stationary policies
                V0(s_,n) = (Para.alpha_1*uAlt(c1_,l1_,Para.psi,Para.sigma)+Para.alpha_2*uAlt(c2_,l2_,Para.psi,Para.sigma))*Para.P(s_,:)'/(1-Para.beta);
                % Store the results for use late as initial guess in the par for loop                
                xInit_0(s_,n,:)=[c1_(s_) c2_(s_) l1_(s_) l2_(s_) x_ R_ x_];
                n=n+1;
                %end
            end
        end
        % Initialize the coeffecients via a routine provided
        % by compEcon library  - funfitxy
        c0(s_,:)=funfitxy(V(s_),squeeze(domain_(s_,:,:)),squeeze(V0(s_,:))' );
    else
        c0(s_,:)=c0(s_-1,:);
        V0(s_,:)=V0(s_-1,:);
        xInit_0(s_,:)=xInit_0(s_-1,:);        
    end
end
domain=vertcat([squeeze(domain_(1,:,:)) ones(length(domain_),1)] ,[squeeze(domain_(1,:,:)) 2*ones(length(domain_),1)]);
c=c0;
save([ Para.datapath 'c1.mat' ] , 'c');

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
    
    


end

