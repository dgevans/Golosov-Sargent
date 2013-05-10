function [PolicyRules, V_new,exitflag,fvec]=CheckGradNAG(xx,RR,ss,c,VV,zInit,Para)
% THIS FUNCTION PERFORMS THE INNER OPTIMIZATION USING THE NAG LIBRARY
% The arguments are explained as follows

%{
xx ,RR,ss : POINT IN THE DOMAIN
c : Coeffs from the current guess
VV: Functional Space
zInit : Initial guess for optimal policies
Para : Parameter strcuts
%}


% THIS ALLOWS US TO USE THE SIMPLE VERSION OF NAG ROOT FINDER.  Technically
% this is no longer necessary but would take to long to change
global V Vcoef R x Par s_ upperFlags lowerFlags

%Get the initial guess for the uconstraint problem. With the simplification
%we need only c1_1,c1_2and c2_1

S = length(Para.P(1,:));
zInit=zInit(1:2*S-1);
Para.alpha=[Para.alpha_1 Para.alpha_2];
%Set global variables
Par=Para;
x=xx;
R=RR;
for s = 1:S
    Vcoef{s}=c(s,:)';
end
V=VV;
s_=ss;
%Set lower and upper limits to the x state variable
xLL=Para.xLL;
xUL=Para.xUL;
n1=Para.n1;
n2=Para.n2;
ctol=Para.ctol;

%% Now solve the unconstraint problem FOC using NAG
% use the last solution
warning('off', 'NAG:warning')
%using nag algorithm find solutions to the FOC
 %[z2, fvec,~,ifail]=c05qb('BelObjectiveUncondGradNAGBGP2Shock',zInit([1,2,4]),'xtol',1e-10);
 %P2shock=[Par.P(1:2,1) sum(Par.P(1:2,2:3),2)];
 %[c12s,c22s,gradc12s,gradc22s] = computeC2_2old(z2(1),z2(2),z2(3),R,s_,P2shock,Par.sigma);
 %zInit=[z2(1) z2(2) z2(2) z2(3) c22s(1,2)];
 %z=zInit;
[z, fvec,~,ifail]=c05qb('BelObjectiveUncondGradNAGBGP',zInit,'xtol',1e-10);
%  if max(abs(z-zInit)) > 1e-7
%      [xx RR]
%      z
%      zInit
%      z-zInit
%  end
%check if code succeeded or failed
       switch ifail
             case {0}
              exitflag=1;
            case {2, 3, 4}
            exitflag=-2;
            z=zInit;
       end


%% GET THE Policy Rules

%Get parameters from Par
psi= Par.psi;
beta =  Par.beta;
P = Par.P;
theta_1 = Par.theta_1;
theta_2 = Par.theta_2;
g = Par.g;
alpha = Par.alpha;
sigma=Par.sigma;
z = z(:)';

c1=z(1:S);
c2_=z(S+1:2*S-1);
  if numel(beta)==1
    beta=ones(2*S-1,S)*beta;
    else
    beta=repmat(beta,2*S-1,1);
    end
%compute components from unconstrained guess
%compute c1 and c2
[c1,c2,gradc1,gradc2] = computeC2_2(c1,c2_,R,s_,P,sigma);
%compute Rprime
[ Rprime,gradRprime ] = computeR( c1,c2,gradc1,gradc2,sigma);
%compute labor supply
[l1 gradl1 l2 gradl2] = computeL(c1,gradc1,c2,gradc2,Rprime,gradRprime,...
                                            theta_1,theta_2,g,n1,n2);
%compute xprime = xprime
[ xprimeMat,gradxprime ] = computeXprime( c1,gradc1,c2,gradc2,Rprime,gradRprime,l1,gradl1,l2,gradl2,...
                                          P,sigma,psi,beta,s_,x);
xprime = xprimeMat(1,:);

% Compute the guess for the multipliers of the constraint problem.
% Lambda_I is multiplier on xprime = xprime (see resFOCBGP_alt.m for
% more detailed description)
dV_x=funeval(Vcoef{1},V(1),[x R],[1 0]);
Lambda_I0=-dV_x;
MultiplierGuess=Lambda_I0 * ones(1,S);

% set flagCons to interior solution
upperFlags = zeros(1,S);
upperFlagsOld= ones(1,S);
lowerFlags = zeros(1,S);
lowerFlagsOld = ones(1,S);

%From solution to unconstrained problem see if upper or lower constraints
%appear to be binding
loopCount=0;
while  (sum(upperFlags~=upperFlagsOld) + sum(lowerFlags~=lowerFlagsOld))>0
    upperFlagsOld = upperFlags;
    lowerFlagsOld = lowerFlags;   
    upperDiff = xprime - xUL;
    upperFlags = ~(upperDiff < 0);
    lowerDiff = xLL - xprime;
    lowerFlags = ~(lowerDiff < 0);
    intFlags = 1- lowerFlags- upperFlags;  
    xdiff = intFlags.*xprime+lowerFlags.*lowerDiff+upperFlags.*upperDiff;
    zInit = [c1(1,:) c2_ xdiff MultiplierGuess];
    %If not in interior
    if sum(intFlags) ~= S
        %% RESOLVE with KKT conditions
        
        warning('off', 'NAG:warning')
        %Find solution to FOCs with extra constraints
      %  SolveKKTWith2Shocks;
       [z, fvec,~,ifail]=c05qb('resFOCBGP_alt',zInit);
        
        z = z(:)';
        
        %Flag if root finding fails
        switch ifail
             case {0}
              exitflag=1;
            case {2, 3, 4}
            exitflag=-2;
            z=zInit;
        end
        
        c1 = z(1:S);
        c2_ = z(S+1:2*S-1);
        
        %compute components from solution
        [c1,c2,gradc1,gradc2] = computeC2_2(c1,c2_,R,s_,P,sigma);
        [ Rprime,gradRprime ] = computeR( c1,c2,gradc1,gradc2,sigma);
        [l1 gradl1 l2 gradl2] = computeL(c1,gradc1,c2,gradc2,Rprime,gradRprime,...
                                            theta_1,theta_2,g,n1,n2);


        
        xprime = intFlags.*z(2*S:3*S-1)+lowerFlags*xLL+upperFlags*xUL;
        loopCount=loopCount+1;
       if loopCount > 3
           disp('Msg:Oscilating xprime')
           disp(xprime)
       end
       if loopCount>5
           break;
       end
    end
    %upperFlagsOld = upperFlags;
    %lowerFlagsOld = lowerFlags;
    
end
%Return policies.
btildprime = xprime./(psi*c2(1,:).^(-sigma));
V_new=-Value3cont(z(1:2*S-1));
%V_new=-Value3cont2Shocks([c1(1,1) c1(1,2) c2(1,1)]);
PolicyRules=[c1(1,:) c2(1,:) l1(1,:) l2(1,:) btildprime Rprime(1,:) xprime];
end
