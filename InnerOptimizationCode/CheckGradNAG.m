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
global V Vcoef R x Par s_ flagCons

%Get the initial guess for the uconstraint problem. With the simplification
%we need only c1_1,c1_2and c2_1
zInit=zInit(1:3);
Para.theta=[Para.theta_1 Para.theta_2];
Para.alpha=[Para.alpha_1 Para.alpha_2];
%Set global variables
Par=Para;
x=xx;
R=RR;
Vcoef{1}=c(1,:)';
Vcoef{2}=c(2,:)';
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
[z, fvec,~,ifail]=c05qb('BelObjectiveUncondGradNAGBGP',zInit,'xtol',1e-10);

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
theta_1 = Par.theta(1);
theta_2 = Par.theta(2);
g = Par.g;
alpha = Par.alpha;
sigma=Par.sigma;
c1_1=z(1);
c1_2=z(2);
c2_1=z(3);

%compute components from unconstrained guess
%compute c1 and c2
[c1,c2,gradc1,gradc2] = computeC2_2(c1_1,c1_2,c2_1,R,s_,P,sigma);
%compute Rprime
[ Rprime,gradRprime ] = computeR( c1,c2,gradc1,gradc2,sigma);
%compute labor supply
[l1 gradl1 l2 gradl2] = computeL(c1,gradc1,c2,gradc2,Rprime,gradRprime,...
                                            theta_1,theta_2,g,n1,n2);
%compute xprime = xprime
[ xprimeMat,gradxprime ] = computeXprime( c1,gradc1,c2,gradc2,Rprime,gradRprime,l1,gradl1,l2,gradl2,...
                                          P,sigma,psi,beta,s_,x);
xprime(1)=xprimeMat(1,1);
xprime(2)=xprimeMat(1,2);

% Compute the guess for the multipliers of the constraint problem.
% Lambda_I is multiplier on xprime = xprime (see resFOCBGP_alt.m for
% more detailed description)
dV_x=[funeval(Vcoef{1},V(1),[x R],[1 0])];
dV_R=[funeval(Vcoef{1},V(1),[x R],[0 1])];
Lambda_I0=-dV_x;
MultiplierGuess=[Lambda_I0 Lambda_I0];
zInit=[c1_1 c1_2 c2_1 xprime(1) xprime(2) MultiplierGuess];

% set flagCons to interior solution
flagCons='ToCheck';
flagConsOld='SolveKKT';

%From solution to unconstrained problem see if upper or lower constraints
%appear to be binding
while  (strcmpi(flagCons,flagConsOld))==0
    flagConsOld=flagCons;
    flagCons='Int';

    % Check the upper limits
    % if upper limit binds for state 1 only
    if xprime(1)> xUL && xprime(2)< xUL
        flagCons='UL_';
        zInit=[c1_1 c1_2 c2_1 (xprime(1)-xUL) xprime(2) MultiplierGuess];

    end
    % if upper limit binds for state 2 only
    if xprime(1) < xUL && xprime(2)>xUL
        flagCons='_UL';
        zInit=[c1_1 c1_2 c2_1 xprime(1)  (xprime(2)-xUL) MultiplierGuess];

    end
    % if upper limit binds for both the states
    if xprime(1)> xUL && xprime(2) > xUL
        flagCons='ULUL';
        zInit=[c1_1 c1_2 c2_1 (xprime(1)- xUL) (xprime(2) - xUL) MultiplierGuess];

    end

    % Check the lower limits
    % if lower limit binds for state 1 only
    if xprime(1)< xLL && xprime(2)> xLL
        flagCons='LL_';
        zInit=[c1_1 c1_2 c2_1 (xLL-xprime(1)) xprime(2) MultiplierGuess];
    end
    % if lower limit binds for state 2 only
    if xprime(1) > xLL && xprime(2) <xLL
        flagCons='_LL';
        zInit=[c1_1 c1_2 c2_1 xprime(1) (xLL-xprime(2)) MultiplierGuess];
    end
    % if lower limit binds for both the states
    if xprime(1) < xLL && xprime(2) <xLL
        flagCons='LLLL';
        zInit=[c1_1 c1_2 c2_1 (xLL-xprime(1)) (xLL-xprime(2)) MultiplierGuess];

    end
    %If not in interior
    if ~(strcmpi(flagCons,'Int'))
        %% RESOLVE with KKT conditions

        warning('off', 'NAG:warning')
        %Find solution to FOCs with extra constraints
        [z, fvec,~,ifail]=c05qb('resFOCBGP_alt',zInit);

        %Flag if root finding fails
        switch ifail
             case {0}
              exitflag=1;
            case {2, 3, 4}
            exitflag=-2;
            z=zInit;
        end

        MuU=zeros(1,2);
        MuL=zeros(1,2);

        c1_1=z(1);
        c1_2=z(2);
        c2_1=z(3);

        %compute components from solution
        [c1,c2,gradc1,gradc2] = computeC2_2(c1_1,c1_2,c2_1,R,s_,P,sigma);
        [ Rprime,gradRprime ] = computeR( c1,c2,gradc1,gradc2,sigma);
        [l1 gradl1 l2 gradl2] = computeL(c1,gradc1,c2,gradc2,Rprime,gradRprime,...
                                                    theta_1,theta_2,g,n1,n2);



        switch flagCons
            case 'LL_'
                % lower limit binds for state 1 only
                MuL(1)=z(4);
                MuL(2)=0;
                xprime(1)=xLL;
                xprime(2)=z(5);

            case '_LL'
                % lower limit binds for state 2 only
                MuL(1)=0;
                MuL(2)=z(5);
                xprime(1)=z(4);
                xprime(2)=xLL;

            case 'LLLL'
                % lower limit binds for both the states
                MuL(1)=z(4);
                MuL(2)=z(5);
                xprime(1)=xLL;
                xprime(2)=xLL;

            case 'UL_'
                % upper limit binds for state 1 only
                MuU(1)=z(4);
                MuU(2)=0;
                xprime(1)=xUL;
                xprime(2)=z(5);

            case '_UL'
                % upper limit binds for state 2 only
                MuU(1)=0;
                MuU(2)=z(5);
                xprime(1)=z(4);
                xprime(2)=xUL;


            case 'ULUL'
                % upper limit binds for both the states
                MuL(1)=z(4);
                MuL(2)=z(5);
                xprime(1)=xUL;
                xprime(2)=xUL;
        end
    end

end
%Return policies.
btildprime = xprime./(psi*c2(1,:).^(-sigma));
V_new=-Value3cont([c1(1,1) c1(1,2) c2(1,1)]);
PolicyRules=[c1(1,1) c1(1,2) c2(1,1) c2(1,2) l1(1,1) l1(1,2) l2(1,1) l2(1,2) btildprime Rprime(1,1) Rprime(1,2) xprime(1) xprime(2)];
end
