% Inputs - xInit, state variables - u2btild,,R,s_  coeff, value
% function, para
function [PolicyRules, V_new,exitflag,fvec]=CheckGradNAG(u2bdiff,RR,s,c,VV,xInit,Para,flagOpt)
global V Vcoef R u2btild Par s_ flagCons

%Get the initial guess for the uconstraint problem. With the simplification
%we need only c1_1,c1_2and c2_1

xInit=xInit(1:3);
Para.theta=[Para.theta_1 Para.theta_2];
Para.alpha=[Para.alpha_1 Para.alpha_2];
Par=Para;
u2btild=u2bdiff;
R=RR;
Vcoef{1}=c(1,:)';
Vcoef{2}=c(2,:)';
V=VV;
s_=s;
u2btildLL=Para.u2btildLL;
u2btildUL=Para.u2btildUL;
n1=Para.n1;
n2=Para.n2;
ctol=Para.ctol;

%% Now solve the unconstraint problem FOC using NAG
% use the last solution
warning('off', 'NAG:warning')
[x, fvec,~,ifail]=c05qb('BelObjectiveUncondGradNAGBGP',xInit);

       switch ifail
             case {0}
              exitflag=1;
            case {2, 3, 4}
            exitflag=-2;
            x=xInit;
        end

if flagOpt==1
    opts = optimset('Algorithm', 'interior-point', 'Display','off','TolX',1e-6);
    xoptguess=x;
    [x, fvec,exitflag]=ktrlink(@(x) -Value3cont(x) ,xoptguess,[],[],[],[],[],[], [],opts);
    exitflag=exitflag+1;
end

%% GET THE Policy Rules
psi= Par.psi;
beta =  Par.beta;
P = Par.P;
theta_1 = Par.theta(1);
theta_2 = Par.theta(2);
g = Par.g;
alpha = Par.alpha;

sigma = 1;
c1_1=x(1);
c1_2=x(2);
c2_1=x(3);

%compute components from unconstrained guess
[c1,c2,gradc1,gradc2] = computeC2_2(c1_1,c1_2,c2_1,R,s_,P,sigma);
[ Rprime,gradRprime ] = computeR( c1,c2,gradc1,gradc2,sigma);
[l1 gradl1 l2 gradl2] = computeL(c1,gradc1,c2,gradc2,Rprime,gradRprime,...
                                            theta_1,theta_2,g,n1,n2);
[ xprime,gradxprime ] = computeXprime( c1,gradc1,c2,gradc2,Rprime,gradRprime,l1,gradl1,l2,gradl2,...
                                          P,sigma,psi,beta,s_,u2btild);
u2btildprime(1)=xprime(1,1);
u2btildprime(2)=xprime(1,2);

% Compute the guess for the multipliers of the constraint problem
dV_x=[funeval(Vcoef{1},V(1),[u2btild R],[1 0])];
dV_R=[funeval(Vcoef{1},V(1),[u2btild R],[0 1])];
Lambda_I0=-dV_x;
MultiplierGuess=[Lambda_I0 Lambda_I0];
xInit=[c1_1 c1_2 c2_1 u2btildprime(1) u2btildprime(2) MultiplierGuess];

% set flagCons to interior solution
flagCons='ToCheck';
flagConsOld='SolveKKT';

while  (strcmpi(flagCons,flagConsOld))==0
    flagConsOld=flagCons;
    flagCons='Int';
    
    % Check the upper limits
    % if upper limit binds for state 1 only
    if u2btildprime(1)> u2btildUL && u2btildprime(2)< u2btildUL
        flagCons='UL_';
        xInit=[c1_1 c1_2 c2_1 (u2btildprime(1)-u2btildUL) u2btildprime(2) MultiplierGuess];
        
    end
    % if upper limit binds for state 2 only
    if u2btildprime(1) < u2btildUL && u2btildprime(2)>u2btildUL
        flagCons='_UL';
        xInit=[c1_1 c1_2 c2_1 u2btildprime(1)  (u2btildprime(2)-u2btildUL) MultiplierGuess];
        
    end
    % if upper limit binds for both the states
    if u2btildprime(1)> u2btildUL && u2btildprime(2) > u2btildUL
        flagCons='ULUL';
        xInit=[c1_1 c1_2 c2_1 (u2btildprime(1)- u2btildUL) (u2btildprime(2) - u2btildUL) MultiplierGuess];
        
    end
    
    % Check the lower limits
    % if lower limit binds for state 1 only
    if u2btildprime(1)< u2btildLL && u2btildprime(2)> u2btildLL
        flagCons='LL_';
        xInit=[c1_1 c1_2 c2_1 (u2btildLL-u2btildprime(1)) u2btildprime(2) MultiplierGuess];
    end
    % if lower limit binds for state 2 only
    if u2btildprime(1) > u2btildLL && u2btildprime(2) <u2btildLL
        flagCons='_LL';
        xInit=[c1_1 c1_2 c2_1 u2btildprime(1) (u2btildLL-u2btildprime(2)) MultiplierGuess];
    end
    % if lower limit binds for both the states
    if u2btildprime(1) < u2btildLL && u2btildprime(2) <u2btildLL
        flagCons='LLLL';
        xInit=[c1_1 c1_2 c2_1 (u2btildLL-u2btildprime(1)) (u2btildLL-u2btildprime(2)) MultiplierGuess];
        
    end
    
    if ~(strcmpi(flagCons,'Int'))
        %% RESOLVE with KKT conditions
        
        warning('off', 'NAG:warning')
        [x, fvec,~,ifail]=c05qb('resFOCBGP_alt',xInit);
        
        switch ifail
             case {0}
              exitflag=1;
            case {2, 3, 4}
            exitflag=-2;
            x=xInit;
        end
        
        MuU=zeros(1,2);
        MuL=zeros(1,2);
        
        c1_1=x(1);
        c1_2=x(2);
        c2_1=x(3);
        
        %compute components from solution
[c1,c2,gradc1,gradc2] = computeC2_2(c1_1,c1_2,c2_1,R,s_,P,sigma);
[ Rprime,gradRprime ] = computeR( c1,c2,gradc1,gradc2,sigma);
[l1 gradl1 l2 gradl2] = computeL(c1,gradc1,c2,gradc2,Rprime,gradRprime,...
                                            theta_1,theta_2,g,n1,n2);


        
        switch flagCons
            case 'LL_'
                % lower limit binds for state 1 only
                MuL(1)=x(4);
                MuL(2)=0;
                u2btildprime(1)=u2btildLL;
                u2btildprime(2)=x(5);
                
            case '_LL'
                % lower limit binds for state 2 only
                MuL(1)=0;
                MuL(2)=x(5);
                u2btildprime(1)=x(4);
                u2btildprime(2)=u2btildLL;
                
            case 'LLLL'
                % lower limit binds for both the states
                MuL(1)=x(4);
                MuL(2)=x(5);
                u2btildprime(1)=u2btildLL;
                u2btildprime(2)=u2btildLL;
                
                
                
            case 'UL_'
                % upper limit binds for state 1 only
                
                MuU(1)=x(4);
                MuU(2)=0;
                u2btildprime(1)=u2btildUL;
                u2btildprime(2)=x(5);
                
                
            case '_UL'
                % upper limit binds for state 2 only
                MuU(1)=0;
                MuU(2)=x(5);
                u2btildprime(1)=x(4);
                u2btildprime(2)=u2btildUL;
                
                
            case 'ULUL'
                
                
                % upper limit binds for both the states
                MuL(1)=x(4);
                MuL(2)=x(5);
                u2btildprime(1)=u2btildUL;
                u2btildprime(2)=u2btildUL;
        end
    end
    
end
btildprime = u2btildprime./(psi*c2(1,:).^(-sigma));
V_new=-Value3cont([c1(1,1) c1(1,2) c2(1,1)]);
PolicyRules=[c1(1,1) c1(1,2) c2(1,1) c2(1,2) l1(1,1) l1(1,2) l2(1,1) l2(1,2) btildprime Rprime(1,1) Rprime(1,2) u2btildprime(1) u2btildprime(2)];
end
