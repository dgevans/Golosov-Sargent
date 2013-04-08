% This program solves recrusive AMSS - optimal policy problem for the
% separable BGP case : u(c,1-n)=psi log (c) + (1-psi) log (1-n)
clc
clear all
close all
opengl software
err=[];
try
    matlabpool('size')
catch err
end

if isempty(err)
   
    if(matlabpool('size') == 0)
        matlabpool open local;
        
    end
    
    
end
err=[];
try
    matlabpool('size')
catch err
end

if isempty(err)
   
    if(matlabpool('size') == 0)
        matlabpool open local;
        
    end
    
    
end

% PARAMETERS
%% Set the technology and preference parameters.
psi=.69; % psi is the relative weight on consmumption in the utility function
beta=.9; % time discount factor
sSize=2; % This is the dimension of the markov shock
g=[.1 .12]; % The vector g is the value of the expenditure shock in each state s
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
ApproxMethod='spli';
OrderOfApprx=20;
orderspli=3;
GridDensity=5;
xGridSize=OrderOfApprx*(GridDensity-1); % size of grid on state variable x
error_tol=1e-6;
NumIter=100;
NumSim=25000;
solveflag=1;
grelax=.9;
Para.flagTransfers=1;
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
compeconpath=[CurrentPath(1:end-length('/AMSS-transfers')) sl 'compecon2011' sl];
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
for s_=1:sSize
n_fb=g+psi*(1-g);
c_fb=psi*(1-g);
uc_fb=1./(1-g);
Euc_fb=sum(pi(s_,:).*uc_fb);
int_fb=uc_fb./(beta*Euc_fb);
x_fb(s_,:)=(((1-psi)/(1)).*((1).*n_fb./(1-n_fb))-psi).*(1-int_fb).^(-1);
end
Para.n_fb=n_fb;
Para.c_fb=c_fb;

xMin=max(x_fb(1,:));
xMax=-xMin;
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
gTrue=Para.g;
g=mean(Para.g)*ones(1,sSize);
n_fb_det=g+psi*(1-g);
c_fb_det=psi*(1-g);
uc_fb_det=1./(1-g);
Euc_fb_det=sum(pi(s_,:).*uc_fb_det);
int_fb_det=uc_fb_det./(beta*Euc_fb_det);
x_fb_det=(((1-psi)/(1)).*((1).*n_fb_det./(1-n_fb_det))-psi).*(1-int_fb_det).^(-1);

for xind=1:xGridSize
    x=xGrid(xind);

for s_=1:sSize
    if ~(x>x_fb_det(s_))
xprime0=ones(1,sSize)*max(x_fb_det(s_));
n(xind,s_,:)=n_fb_det;
u(xind,s_,:)=psi*log(squeeze(n(xind,s_,:))'-g)+(1-psi)*log(1-squeeze(n(xind,s_,:))');
Eu(xind,s_)=sum(pi(s_,:).*squeeze(u(xind,s_,:))');
    else
xprime0=[x x];
n(xind,s_,:)=nDet(x);
u(xind,s_,:)=psi*log(squeeze(n(xind,s_,:))'-g)+(1-psi)*log(1-squeeze(n(xind,s_,:))');
Eu(xind,s_)=sum(pi(s_,:).*squeeze(u(xind,s_,:))');
    
    end
n0guess=squeeze(n(xind,s_,:))';
end
end
V0=(inv(A)*Eu')';
for s=1:sSize
coeff0(s,:)=funfitxy(V(s),xGrid,V0(:,s));
end

Err0=sqrt(sum((funeval(coeff0(1,:)',V(1),xGrid)-V0(:,1)).^2));
plot(xGrid,V0)
g=gTrue;

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
parfor xind=1:xGridSize
    x=xGrid(xind);
    s_=1;
%for s_=1:sSize
    [nguess,xprimeguess,c,~] =solveInnerOpt(x,s_,coeff,V,Para);
  [n(xind,1,:),xprime(xind,1,:),c,VNew(xind,1),exitflag(xind,1,:)] =solveInnerOptUsingConsOptFminCon(x,s_,coeff,V,Para,[nguess,xprimeguess]);
 
%end
end
  n(:,2,:)=n(:,1,:);
  xprime(:,2,:)=xprime(:,1,:);
  VNew(:,2)=VNew(xind,1);
  exitflag(:,2,:)=exitflag(:,1,:);
 
%% Update the coeff
coeffold=coeff;
for s=1:sSize
coeffNew(s,:)=funfitxy(V(s),xGrid(logical((exitflag(:,s_,:)==solveflag))),VNew(logical((exitflag(:,s_,:)==solveflag))));
%coeffNew(s,:)=funfitxy(V(s),xGrid,VNew);
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
 
save('~/Golosov-Sargent/Data/temp/AMSS_Solution_tr.mat','coeff','V','Para','n','xprime','Err','coeffN','N','coeffXPrime','XPrime')



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


NumSim=25000
rhist0=rand(NumSim,1);
s0=1
b_=0;
[SimDataPol]=runSimulationUsingPolicyRules(b_,s0,Para,coeffN,N,coeff,V,rhist0,NumSim)
%[SimDataOpt]=runSimulation(b_,s0,Para,coeff,V,rhist0,NumSim);
figure()
plot(SimDataPol.xHist)
print(gcf,'-dpng',[ plotpath 'FigFigAMSSSimulation.png'])

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
