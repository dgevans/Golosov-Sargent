function [c,x_state,PolicyRulesStore]=InitializeWithDeterministicValues(V,Para,flagComputeInitCoeff)

gTrue=Para.g;
Para.g=mean(gTrue)*ones(2,1);
n1=Para.n1;
    n2=Para.n2;
    alpha_1=Para.alpha_1;
    alpha_2=Para.alpha_2;
    g=Para.g(1);
    theta_1=Para.theta_1;
    theta_2=Para.theta_2;
    psi=Para.psi;
    sigma=Para.sigma;
gTrue=Para.g;
Para.g=mean(gTrue)*ones(2,1);
   
for s_=1:Para.sSize
    n=1;
    if s_==1
        
        
        for u2btildctr=1:Para.u2btildGridSize
            for Rctr=1:Para.RGridSize
                
                u2btild_=Para.u2bdiffGrid(u2btildctr);
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
            c2_=R_^(-1/Para.sigma)*c1_;
                    % Solve for l1 , l2 using the resource constraint and wage equation
                TotalResources=(c1_*n1+c2_*n2+g);
DenL2=theta_2*R_*n1+theta_2*n2;
l2_=(TotalResources-theta_1*n1+ theta_2*n1*R_)/(DenL2);
if theta_2==0
    l1_=TotalResources/(n1*theta_1);
l2_=0;
else
l1_= 1-(1-l2_)*theta_2/theta_1*R_;
end

                u2btildPrime_=u2btild_;
                V0(s_,n)=(Para.alpha_1*uAlt(c1_,l1_,Para.psi,Para.sigma)+Para.alpha_2*uAlt(c2_,l2_,Para.psi,Para.sigma))/(1-Para.beta);
                if strcmpi(flagComputeInitCoeff,'no')
                    V0(s_,n)=funeval(cInit(s_,:)',VInit(s_),[u2btild_,R_]);
                    PolicyRulesStore(n,:)=GetInitialApproxPolicy([u2btild_,R_,s_],InitData.x_state,InitData.PolicyRulesStore);
                end
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
        if strcmpi(flagComputeInitCoeff,'no')
            PolicyRulesStore(GridSize/2+1:GridSize,:)=PolicyRulesStore(1:GridSize/2,:);
        end
        
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

if strcmpi(flagComputeInitCoeff,'yes')
    PolicyRulesStore=vertcat(PolicyRulesStore1,PolicyRulesStore2);
end