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

%% Set the technology and preference parameters.
psi=.7; % psi is the relative weight on consmumption in the utility function
beta=.9; % time discount factor
g=[.1 .15]; % The vector g is the value of the expenditure shock in each state s
sSize=length(g); % This is the dimension of the markov shock
pi=ones(sSize,sSize)/sSize;
%pi=repmat([.7 .1 .2],sSize,1);
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
OrderOfApprx=25;
orderspli=3;
GridDensity=2;
xGridSize=OrderOfApprx*(GridDensity-1); % size of grid on state variable x
error_tol=1e-7;
NumIter=200;
NumSim=25000;
solveflag=0;
Para.solveflag=solveflag;
grelax=.9;
Para.flagTransfers=1;
no_transfer_max_iter=1;
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
Para.x_fb=min(x_fb(1,:));
xMin=min(x_fb(1,:));
xMax=-xMin*.9;
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
g=sum(Para.g.*pi(1,:))*ones(1,sSize);
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
xprime0=x*ones(1,sSize);
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

% Initialize with the no-transfers case
%AMSSNoTransfers=load('~/Golosov-Sargent/Data/temp/AMSS_Solution_no_tr_2shocks.mat');

%for s=1:sSize
%coeff0(s,:)=funfitxy(V(s),xGrid,funeval(AMSSNoTransfers.coeff(s,:)',AMSSNoTransfers.V(s),xGrid ));
%end


%% VALUE FUNCTION ITERATION
options = optimset('Display','off','TolX', 1e-10,'TolCon',1e-10,'TolFun',1e-10);
Para.options=options;
n=ones(xGridSize,sSize,sSize)*.5;
coeff=coeff0;
iter=1;
error=1;
while( error>error_tol && iter<NumIter)
tic
s_=1;
parfor xind=1:xGridSize
    x=xGrid(xind);
    
    [n(xind,s_,:),xprime(xind,s_,:),c,VNew(xind,s_),exitflag(xind,s_,:)] =solveInnerOpt(x,s_,coeff,V,Para);
    
    if iter>no_transfer_max_iter
   [n(xind,s_,:),xprime(xind,s_,:),c,VNew(xind,s_),exitflag(xind,s_,:)] =...
        solveInnerOptUsingConsOptFminCon(x,s_,coeff,V,Para,[squeeze(n(xind,s_,:))',squeeze(xprime(xind,s_,:))']);
    end
end

for s_=2:sSize
  n(:,s_,:)=n(:,1,:);
  xprime(:,s_,:)=xprime(:,1,:);
  VNew(:,s_)=VNew(:,1);
  exitflag(:,s_,:)=exitflag(:,1,:);
end
%% Update the coeff
coeffold=coeff;
for s=1:sSize
coeffNew(s,:)=funfitxy(V(s),xGrid(logical((exitflag(:,s_,:)==solveflag))),VNew(logical((exitflag(:,s_,:)==solveflag))));
%coeffNew(s,:)=funfitxy(V(s),xGrid,VNew);
[ coeffNew(s,:)] = SmoothPasting(V(s),VNew(:,s),xGrid,x_fb(1,1),[]);
ValueDiff(s,:)=funeval(coeffold(s,:)',V(s),xGrid)-funeval(coeffNew(s,:)',V(s),xGrid);
end
coeff=coeffold*(1-grelax)+coeffNew*grelax;
Err(iter,:)=((sum(ValueDiff.^2,2)).^(.5))';
error=Err(iter)
iter=iter+1
toc
save('~/Golosov-Sargent/Data/temp/AMSS_Solution_tr_2shocks_bgp.mat','coeff','V','Para','n','xprime','Err','xGrid')
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
 
save('~/Golosov-Sargent/Data/temp/AMSS_Solution_tr_2shocks_bgp.mat','coeff','V','Para','n','xprime','Err','coeffN','N','coeffXPrime','XPrime','xGrid')


break;



    
