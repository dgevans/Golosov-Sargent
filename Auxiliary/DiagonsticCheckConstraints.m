% Inputs - xInit, state variables - u2btild,,R,s_  coeff, value
% function, para
function [Diagnostic1]=DiagonsticCheckConstraints(u2bdiff,RR,s,c,VV,xInit,Para)
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
[x, fvec,exitflag]=c05nb('BelObjectiveUncondGradNAGBGP',xInit,'xtol',1e-10)

h=.00001
Der(1)= (Value3cont([x(1)+h,x(2),x(3)])-Value3cont([x(1)-h,x(2),x(3)]))/(2*h);
Der(2)= (Value3cont([x(1),x(2)+h,x(3)])-Value3cont([x(1),x(2)-h,x(3)]))/(2*h);
Der(3)= (Value3cont([x(1),x(2),x(3)+h])-Value3cont([x(1),x(2),x(3)-h]))/(2*h);
Der


%opts = optimset('Algorithm', 'interior-point', 'Display','iter','TolX',1e-6);    
xoptguess=x;
[x, fvec,exitflag]=ktrlink(@(x) -Value3cont(x) ,xoptguess,[],[],[],[],[],[], [],[]);
x-xoptguess

figure()
 subplot(1,3,1)
 fplot (@(x1) Value3cont([x1 x(2) x(3)]),[x(1)*.8 x(1)*1.25] )
 subplot(1,3,2)
 fplot (@(x2) Value3cont([x(1) x2 x(3)]),[x(2)*.8 x(2)*1.25] )
 subplot(1,3,3)
 fplot (@(x3) Value3cont([x(1) x(2) x3]),[x(3)*.8 x(3)*1.25] )


% use a derivative free optimizer



if exitflag==4
    exitflag=-2;
else
    exitflag=1;
end
psi= Par.psi;
beta =  Par.beta;
P = Par.P;
theta_1 = Par.theta(1);
theta_2 = Par.theta(2);
g = Par.g;
alpha = Par.alpha;

sigma = 1;
    frac = (R*P(s_,1)*x(1)^(-sigma)+R*P(s_,2)*x(2)^(-sigma)-P(s_,1)*x(3)^(-sigma))...
        /( P(s_,2) );
     c1_1=x(1);
     c1_2=x(2);
     c2_1=x(3);
     
     
     
    %compute components from unconstrained guess
frac = (R*P(s_,1)*c1_1^(-sigma)+R*P(s_,2)*c1_2^(-sigma)-P(s_,1)*c2_1^(-sigma))...
        /( P(s_,2) );
    c2_2 = frac^(-1/sigma);
    
     %Compute l1 form formula
    l1_1den = n1*theta_1+n2*c2_1*theta_1/c1_1;
    l1_1num = (n1*c1_1+n2*c2_1+g(1) + n2*(c2_1*theta_1-c1_1*theta_2)/c1_1);
    l1(1) = l1_1num/l1_1den;
    l1_2den = n1*theta_1+n2*c2_2*theta_1/c1_2;
    l1_2num = (n1*c1_2+n2*c2_2+g(2) + n2*(c2_2*theta_1-c1_2*theta_2)/c1_2);
    l1(2) = l1_2num/l1_2den;
      %compute l2 from formula
    l2_1den = n2*theta_2+n1*c1_1*theta_2/c2_1;
    l2_1num = n1*c1_1+n2*c2_1+g(1)+n1*(c1_1*theta_2-c2_1*theta_1)/c2_1;
    l2(1) = l2_1num/l2_1den;
    l2_2den = n2*theta_2+n1*c1_2*theta_2/c2_2;
    l2_2num = n1*c1_2+n2*c2_2+g(2)+n1*(c1_2*theta_2-c2_2*theta_1)/c2_2;
    l2(2) = l2_2num/l2_2den;
    %get expected value of marginal utility of agent 2
    Eu2 = P(s_,1)*c2_1^(-1)+P(s_,2)*c2_2^(-1);
    
    %compute btildeprime from formula
    btildprime(1) = u2btild/(beta*psi*Eu2)...
        +c1_1-c2_1-(1-psi)*c1_1*l1(1)/(psi*(1-l1(1)))+(1-psi)*c2_1*l2(1)/(psi*(1-l2(1)));

    %Compute btildprime(2) from formula
    btildprime(2) = u2btild/(beta*psi*Eu2)...
        +c1_2-c2_2-(1-psi)*c1_2*l1(2)/(psi*(1-l1(2)))+(1-psi)*c2_2*l2(2)/(psi*(1-l2(2)));

u2btildprime=psi*[c2_1^(-1) c2_2^(-1)].*btildprime;


X(1,:) = [psi*c2_1^(-1)*btildprime(1),c2_1^(-1)/c1_1^(-1)];%state next period
X(2,:) = [psi*c2_2^(-1)*btildprime(2),c2_2^(-1)/c1_2^(-1)];%state next period

% Compute the guess for the multipliers of the constraint problem
dV_x=[funeval(Vcoef{1},V(1),[u2btild R],[1 0])];
dV_R=[funeval(Vcoef{1},V(1),[u2btild R],[0 1])];
Lambda_I0=-dV_x;
MultiplierGuess=[Lambda_I0 Lambda_I0];

%dV_x=[funeval(Vcoef{1},V(1),[u2btild R],[1 0])];
%dV_R=[funeval(Vcoef{1},V(1),[u2btild R],[0 1])];
%Lambda_I0=-dV_x*u2btild*beta;
%Lambda_B0=-dV_R;
%Lambda_R0 = -P(s_,:).*psi.*[c1_1^(-1) c1_2^(-1)];
%Lambda_W0 = -P(s_,:).*[-(1-psi)/(1-l1(1)) + theta_1*psi*c2_1^(-1) -(1-psi)/(1-l1(2)) + theta_1*psi*c2_2^(-1)];
%MultiplierGuess=[Lambda_I0 Lambda_I0 Lambda_B0 Lambda_R0 Lambda_W0];
%xInit=[c1_1 c1_2 c2_1 c2_2 l1 l2 u2btildprime(1) u2btildprime(2) MultiplierGuess];
xInit=[c1_1 c1_2 c2_1 u2btildprime(1) u2btildprime(2) MultiplierGuess];

% set flagCons to interior solution
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
 flagCons='Int';
if ~(strcmpi(flagCons,'Int'))
    
    warning('off', 'NAG:warning')
[x, fvec,exitflag]=c05nb('resFOCBGP_alt',xInit,'xtol',1e-8);
%[x, fvec,exitflag]=c05nb('resFOC',xInit);
if exitflag==4
    exitflag=-2;
    %x=xInit;
else
    exitflag=1;
end
    % solve for the constrainted problem
  %  options = optimset('Display','off','TolFun',ctol,'FunValCheck','off','TolX',ctol,'MaxFunEvals', 100*length(xInit),'MaxTime',100);
   % [x fval exitflag] =fsolve(@(x) resFOC(x,u2btild,R,c,s_,V,flagCons,Para),xInit,options)  ;
    MuU=zeros(1,2);
MuL=zeros(1,2);
   frac = (R*P(s_,1)*x(1)^(-sigma)+R*P(s_,2)*x(2)^(-sigma)-P(s_,1)*x(3)^(-sigma))...
        /( P(s_,2) );
     c1_1=x(1);
     c1_2=x(2);
     c2_1=x(3);
     
     
     
    %compute components from unconstrained guess
frac = (R*P(s_,1)*c1_1^(-sigma)+R*P(s_,2)*c1_2^(-sigma)-P(s_,1)*c2_1^(-sigma))...
        /( P(s_,2) );
    c2_2 = frac^(-1/sigma);
    
     %Compute l1 form formula
    l1_1den = n1*theta_1+n2*c2_1*theta_1/c1_1;
    l1_1num = (n1*c1_1+n2*c2_1+g(1) + n2*(c2_1*theta_1-c1_1*theta_2)/c1_1);
    l1(1) = l1_1num/l1_1den;
    l1_2den = n1*theta_1+n2*c2_2*theta_1/c1_2;
    l1_2num = (n1*c1_2+n2*c2_2+g(2) + n2*(c2_2*theta_1-c1_2*theta_2)/c1_2);
    l1(2) = l1_2num/l1_2den;
      %compute l2 from formula
    l2_1den = n2*theta_2+n1*c1_1*theta_2/c2_1;
    l2_1num = n1*c1_1+n2*c2_1+g(1)+n1*(c1_1*theta_2-c2_1*theta_1)/c2_1;
    l2(1) = l2_1num/l2_1den;
    l2_2den = n2*theta_2+n1*c1_2*theta_2/c2_2;
    l2_2num = n1*c1_2+n2*c2_2+g(2)+n1*(c1_2*theta_2-c2_2*theta_1)/c2_2;
    l2(2) = l2_2num/l2_2den;
    %get expected value of marginal utility of agent 2
    Eu2 = P(s_,1)*c2_1^(-1)+P(s_,2)*c2_2^(-1);
    
    %compute btildeprime from formula
    btildprime(1) = u2btild/(beta*psi*Eu2)...
        +c1_1-c2_1-(1-psi)*c1_1*l1(1)/(psi*(1-l1(1)))+(1-psi)*c2_1*l2(1)/(psi*(1-l2(1)));

    %Compute btildprime(2) from formula
    btildprime(2) = u2btild/(beta*psi*Eu2)...
        +c1_2-c2_2-(1-psi)*c1_2*l1(2)/(psi*(1-l1(2)))+(1-psi)*c2_2*l2(2)/(psi*(1-l2(2)));


    
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
        
       otherwise
       MuL(1)=0;
       MuL(2)=0;
       u2btildprime(1)=x(4);
       u2btildprime(2)=x(5);     
       
        
    end
   
    X(1,:) = [u2btildprime(1),c2_1^(-1)/c1_1^(-1)];%state next period
    X(2,:) = [u2btildprime(2),c2_2^(-1)/c1_2^(-1)];%state next period
    Lambda=x(11:end);
   
%end





%compute objective
if ~isreal(X)
    X=real(X);
end
Vobj = P(s_,1)*(alpha(1)*uBGP(c1_1,l1(1),psi)+alpha(2)*uBGP(c2_1,l2(1),psi)...
    +beta*funeval(Vcoef{1},V(1),X(1,:)));

Vobj = Vobj + P(s_,2)*(alpha(1)*uBGP(c1_2,l1(2),psi)+alpha(2)*uBGP(c2_2,l2(2),psi)...
    +beta*funeval(Vcoef{2},V(2),X(2,:)));

V_new=Vobj;
PolicyRules=[c1_1 c1_2 c2_1 c2_2 l1(1) l1(2) l2(1) l2(2) btildprime c2_1^(-1)/c1_1^(-1) c2_2^(-1)/c1_2^(-1) u2btildprime(1) u2btildprime(2)];
 [ EqCons,EqConsJac ] =NonLinearConstraints([c1_1 c1_2 c2_1 c2_2 l1(1) l1(2) l2(1) l2(2) u2btildprime(1) u2btildprime(2)]);
 Diagnostic1=max(EqCons);
 if Diagnostic1 >1e-6
     disp('Error')
     disp(Diagnostic1)
     disp(u2bdiff)
     disp(RR)
     disp(flagCons)
     exitflag
 end
end

