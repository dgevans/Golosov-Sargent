% This program solves recrusive AMSS - optimal policy problem for the
% separable BGP case : u(c,1-n)=psi log (c) + (1-psi) log (1-n)
clc
clear all
close all
opengl software
% PARAMETERS
%% Set the technology and preference parameters.
psi=.69; % psi is the relative weight on consmumption in the utility function
beta=.9; % time discount factor
sSize=2; % This is the dimension of the markov shock
g=[.1 .15]; % The vector g is the value of the expenditure shock in each state s
pi=[.5 .5;
    .5 .5];  % The probability transition matrix for the markov shock
A=eye(sSize);
for s=1:sSize
    for sprime=1:sSize
        if ~(s==sprime)
        A(s,sprime)=-beta*pi(s,sprime);
        else
                   A(s,sprime)=1-beta*pi(s,sprime);
        end 
    end
end

%% Set other parameters
Para.psi=psi;
Para.g=g;
Para.pi=pi;
Para.beta=beta;
Para.sSize=sSize;
Para.invA=inv(A);
ApproxMethod='cheb';
OrderOfApprx=20;
orderspli=3;
GridDensity=4;
xGridSize=OrderOfApprx*(GridDensity-1); % size of grid on state variable x
error_tol=1e-5;
NumIter=100;
NumSim=200;
solveflag=0;
grelax=.9;

% PATHS
% Compecon

if strcmp(computer,'PCWIN')
    sl='\';
    coresize=2;
else
    sl='/';
    coresize=8;
end
CurrentPath=pwd;
compeconpath=[CurrentPath(1:end-length('/AMSS')) sl 'compecon2011' sl];
texpath= [pwd sl 'Tex' sl] ;
plotpath= [pwd sl 'Graphs' sl] ;
datapath=[pwd sl 'Data' sl] ;
mkdir(texpath)
mkdir(plotpath)
mkdir(datapath)
addpath(genpath(compeconpath))
Para.plotpath=plotpath;
%%
% GRID AND THE FUNCTIONAL SPACE
% We now compute the grid for x. 
negtaxrev=@(n) -beta*(n*(1-((1-psi)/(psi))*((n-min(g))/(1-n)))-min(g))/(1-beta);
[~,xMax]=fminbnd(negtaxrev,min(g),1);
xMax=-xMax;
xMin=-(beta/(1-beta))*(max(g)/(1-max(g)));
xMax=2;
xMin=-2;
% Define the functional space
for s=1:sSize
V(s) = fundefn(ApproxMethod,OrderOfApprx ,xMin,xMax,orderspli);
end
xGrid=linspace(xMin,xMax,xGridSize)';

xGridSize=length(xGrid);


Para.xMax=xMax;
Para.xMin=xMin;
Para.xGrid=xGrid;
%% Initialize the value function
nDet=@(x) (psi+x*((beta-1)/beta))/(1+x*((beta-1)/beta));
Para.nDet=nDet;
for xind=1:xGridSize
    x=xGrid(xind);

for s_=1:sSize
xprime0=[x x];
n(xind,s_,:)=nDet(x);
u(xind,s_,:)=psi*log(squeeze(n(xind,s_,:))'-g)+(1-psi)*log(1-squeeze(n(xind,s_,:))');
Eu(xind,s_)=sum(pi(s_,:).*squeeze(u(xind,s_,:))');
n0guess=squeeze(n(xind,s_,:))';
end
end
V0=(inv(A)*Eu')';
for s=1:sSize
coeff0(s,:)=funfitxy(V(s),xGrid,V0(:,s));
end

Err0=sqrt(sum((funeval(coeff0(1,:)',V(1),xGrid)-V0(:,1)).^2));


%% VALUE FUNCTION ITERATION
options = optimset('Display','off');
Para.options=options;
nguess=[nDet(xGrid(1)) nDet(xGrid(2))];
n=ones(xGridSize,sSize,sSize)*.5;
coeff=coeff0;
iter=1;
error=1;
while( error>error_tol && iter<NumIter)
tic
for xind=1:xGridSize
    x=xGrid(xind);
for s_=1:sSize
  [n(xind,s_,:),xprime(xind,s_,:),c,exitflag(xind,s_,:)] =solveInnerOpt(x,s_,coeff,V,Para);

%% Compute the fitted points 
for s=1:sSize
Vstar(s)=funeval(coeff(s,:)',V(s),xprime(xind,s_,s));
end
u(xind,s_,:)=psi*log(squeeze(n(xind,s_,:))'-g)+(1-psi)*log(1-squeeze(n(xind,s_,:))')+beta*Vstar;
VNew(xind,s_)=sum(pi(s_,:).*squeeze(u(xind,s_,:))');
nguess=squeeze(n(xind,s_,:))';
end
end
%% Update the coeff
coeffold=coeff;
for s=1:sSize
coeffNew(s,:)=funfitxy(V(s),xGrid(logical((exitflag(:,s_,:)==solveflag))),VNew(logical((exitflag(:,s_,:)==solveflag))));
ValueDiff(s,:)=funeval(coeffold(s,:)',V(s),xGrid)-funeval(coeffNew(s,:)',V(s),xGrid);
end
coeff=coeffold*(1-grelax)+coeffNew*grelax;
Err(iter,:)=((sum(ValueDiff.^2,2)).^(.5))';
error=Err(iter)
iter=iter+1
toc
end


% approximate policy rules
for s=1:sSize
N(s) = fundefn(ApproxMethod,OrderOfApprx ,xMin,xMax,orderspli);
XPrime(s)=fundefn(ApproxMethod,OrderOfApprx ,xMin,xMax,orderspli);
coeffN(s,:)=funfitxy(N(s),xGrid,n(:,1,s));
coeffXPrime(s,:)=funfitxy(XPrime(s),xGrid,xprime(:,1,s));
end




 
% %% Time 0 problem
 for s=1:sSize
     s0=s;
 for xind=1:xGridSize
     
 b_=xGrid(xind)*(1-g(s0));
 n0guess=nDet((1/1-g(s0))*b_);
     get_root_labor0_nag= @(num,n0,user,iflag) getResLaborTime0(num,b_,user,s0,n0,coeff,V,Para,iflag);
 [n0(xind,s0),~,exitflag0]=c05qb(get_root_labor0_nag,n0guess,'xtol',1e-10);
 ExitFlag0(xind,s0)=exitflag0;
 c0(xind,s0)=n0(xind,s0)-g(s);
 xprime0(xind,s0)=(n0(xind,s0)/(1-n0(xind,s0)))*(1-psi)+ (psi/c0(xind,s0))*b_-psi;
 n0guess=n0(xind,s0);
 end
 end
 
save('~/Golosov-Sargent/Data/temp/AMSS_Solution.mat','coeff','V','Para','n','xprime','Err','coeffN','N','coeffXPrime','XPrime')



% Figure 1 : Value function
figure()
s_=1;
plot(xGrid,funeval(coeff(s_,:)',V(s_),xGrid),':k','LineWidth',2)
hold on
s_=2;
plot(xGrid,funeval(coeff(s_,:)',V(s_),xGrid),'k','LineWidth',2)
xlabel('x')
ylabel('V(x)')
print(gcf,'-dpng',[ plotpath 'FigAMSSValueFunction1.png'])

% Figure 2a : Policy Rule -n
figure()
s_=1;
subplot(1,2,1)
plot(xGrid,squeeze(n(:,s_,1)),':k','LineWidth',2)
hold on
plot(xGrid,squeeze(n(:,s_,2)),'k','LineWidth',2)
xlabel('x')
ylabel('n(x)')
title('$s_{-1}=1$','Interpreter','Latex')

subplot(1,2,2)
s_=2;
plot(xGrid,squeeze(n(:,s_,1)),':k','LineWidth',2)
hold on
plot(xGrid,squeeze(n(:,s_,2)),'k','LineWidth',2)
xlabel('x')
ylabel('n(x)')
title('$s_{-1}=2$','Interpreter','Latex')
print(gcf,'-dpng',[ plotpath 'FigAMSSLaborPolicy1.png'])

% Figure 2a : Policy Rule - taxes
% tau=1-ul/uc = 1- (1-psi)/()
tax=@(n,s) 1-((1-psi)/(psi)).*(n-g(s))./(1-n);
figure()
s_=1;
subplot(1,2,1)
plot(xGrid,tax(squeeze(n(:,s_,1)),1),':k','LineWidth',2)
hold on
plot(xGrid,tax(squeeze(n(:,s_,2)),2),'k','LineWidth',2)
xlabel('x')
ylabel('\tau(x)')
title('$s_{-1}=1$','Interpreter','Latex')

subplot(1,2,2)
s_=2;
plot(xGrid,tax(squeeze(n(:,s_,1)),1),':k','LineWidth',2)
hold on
plot(xGrid,tax(squeeze(n(:,s_,2)),2),'k','LineWidth',2)
xlabel('x')
ylabel('\tau(x)')
title('$s_{-1}=1$','Interpreter','Latex')
print(gcf,'-dpng',[ plotpath 'FigAMSSTaxPolicy.png'])


% Figure 2b : Policy Rule -x'

 figure()
 s_=1;
 subplot(1,2,1)
 plot(xGrid,squeeze(xprime(:,s_,1))-xGrid,':k','LineWidth',2)
 hold on
 plot(xGrid,squeeze(xprime(:,s_,2))-xGrid,'k','LineWidth',2)
 hold on

 s_=2;
 subplot(1,2,2)
 plot(xGrid,squeeze(xprime(:,s_,1))-xGrid,':k','LineWidth',2)
 hold on
 plot(xGrid,squeeze(xprime(:,s_,2))-xGrid,'k','LineWidth',2)
 hold on
 print(gcf,'-dpng',[ plotpath 'FigAMSSxPolicy1.png'])
 
 
 %% Diagnostics
 figure()
 % This figure plots the L2 norm error for the value function convergence
 % across iterations
 plot(Err)
 hold on
 plot((1:NumIter),repmat([1e-3 -1e-3],NumIter,1),':k','LineWidth',2)
 axis([1 NumIter -1e-4 1e-2])
 
% error in policy rule approximations 
figure()
s=1;
plot(xGrid,funeval(coeffN(s,:)',N(s),xGrid)-n(:,1,s),'k','LineWidth',2)
hold on
s=1;
plot(xGrid,funeval(coeffN(s,:)',N(s),xGrid)-n(:,1,s),':k','LineWidth',2)
xlabel('x')
ylabel('Err(x)')
print(gcf,'-dpng',[ plotpath 'FigAMSSErrorPolicyRules_n.png'])


figure()
s=1;
plot(xGrid,funeval(coeffXPrime(s,:)',XPrime(s),xGrid)-xprime(:,1,s),'k','LineWidth',2)
hold on
s=1;
plot(xGrid,funeval(coeffXPrime(s,:)',XPrime(s),xGrid)-xprime(:,1,s),':k','LineWidth',2)
xlabel('x')
ylabel('Err(x)')
print(gcf,'-dpng',[ plotpath 'FigAMSSErrorPolicyRules_xprime.png'])


NumSim=10000
rhist0=rand(NumSim,1);
s0=1
b_=-1;
[SimDataPol]=runSimulationUsingPolicyRules(b_,s0,Para,coeffN,N,coeff,V,rhist0,NumSim)
%[SimDataOpt]=runSimulation(b_,s0,Para,coeff,V,rhist0,NumSim);
figure()
plot(SimDataPol.xHist)
%max(SimDataPol.xHist./SimDataOpt.xHist)-1
%% TO DO
% Check the bounds on xprime
% generalize for ces, ghh and other preferences
% extra code for grid
%  x=mean(xGrid(end));
%  s_=1;
%  nguess=[nDet(x) nDet(x)];
%  delta=.001
%  [DerL]=.5*checkGradients(x,s_,nguess,coeff,V,Para,delta)+.5*checkGradients(x,s_,nguess,coeff,V,Para,-delta)
%  [resINQ,resEQ]=getResLaborKnitro(x,s_,nguess,coeff,V,Para)
% xGridDefault=funnode(V(1));
% xGrid=ones(xGridSize,1);
% xGrid(1)=xGridDefault(1);
% for xind=2:xGridSize
%     if mod(xind,GridDensity)==0
%         if xind==GridDensity
%     xGrid(xind-GridDensity+1:xind)=linspace(xGridDefault(xind/GridDensity),xGridDefault(xind/GridDensity+1),GridDensity);
%         else
%     xGrid(xind-GridDensity:xind-1)=linspace(xGridDefault(xind/GridDensity),xGridDefault(xind/GridDensity+1),GridDensity);
%             
%     end
%     end
% end
%xGrid=linspace(xMin,xMax*.9,2*xGridSize/3)';
%xGrid=[xGrid ;linspace(xMax*.91,xMax,xGridSize/3)'];
