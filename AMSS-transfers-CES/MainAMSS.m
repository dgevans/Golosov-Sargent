% This program solves recrusive AMSS - optimal policy problem for the
% separable BGP case : u(c,1-n)=psi log (c) + (1-psi) log (1-n)
function MainAMSS(no_transfer_max_iter,storefilename,g,psi,sigma,gamma,beta)
if(matlabpool('size') == 0)
matlabpool open local;
end 

%% Set the technology and preference parameters.
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

%Populate the utility function and the required derivatives
if ~(sigma ==1)
util=@(n) (n-g).^(1-sigma)./(1-sigma)-n.^(1+gamma)./(1+gamma); % u(n-g,n)
else
util=@(n) log(n-g)-n.^(1+gamma)./(1+gamma); % u(n-g,n)
end
der_u_c=@(n) (n-g).^(-sigma); % d u /d c (as a function of n)
der_u_n=@(n) n.^(gamma); % d u / d n
der_u_cn=@(n) -sigma*(n-g).^(-sigma-1); % d u /dc^2
der_u_nn=@(n) gamma*n.^(gamma-1); % d u / dn^2
res_imp_det =@(n,x) (n-sum(g.*pi(1,:)))^(-sigma)*(n-sum(g.*pi(1,:)))+x-n^(1+gamma)-x/beta;
 

%%BGP
%util=@(n) psi*log(n-g)+(1-psi)*log(1-n); % u(n-g,n)
%der_u_c=@(n) psi*(n-g).^(-1); % d u /d c (as a function of n)
%der_u_n=@(n) (1-psi)*(1-n).^(-1); % d u / d n
%der_u_cn=@(n) -1*(n-g).^(-2); % d u /dc^2
%der_u_nn=@(n) (1-psi)./(1-n).^2; % d u / dn^2
%res_imp_det =@(n,x) psi*(n-sum(g.*pi(1,:)))^(-1)*(n-sum(g.*pi(1,:)))+x- (1-psi)*(n/(1-n)) -x/beta;


Para.der_u_c=der_u_c;
Para.der_u_n=der_u_n;
Para.util=util;
Para.der_u_cn=der_u_cn;
Para.der_u_nn=der_u_nn;
%% Set other parameters
nDet=@(x) (psi+x*((beta-1)/beta))/(1+x*((beta-1)/beta));
Para.nDet=nDet;
Para.psi=psi;
Para.sigma=sigma;
Para.gamma=gamma;
Para.g=g;
Para.pi=pi;
Para.beta=beta;
Para.sSize=sSize;
Para.invA=inv(A);
Para.der_u_c=der_u_c;
Para.der_u_n=der_u_n;
Para.util=util;
Para.u_cn=der_u_cn;
Para.u_nn=der_u_nn;
ApproxMethod='spli';
OrderOfApprx=20;
orderspli=3;
GridDensity=3;
xGridSize=OrderOfApprx*(GridDensity-1); % size of grid on state variable x
error_tol=1e-7;
NumIter=100;
solveflag=0;
Para.solveflag=solveflag;
grelax=.95;
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
compeconpath=[CurrentPath(1:end-length('/AMSS-transfers-CES')) sl 'compecon2011' sl];
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
[n_fb(s_,:),c_fb(s_,:),x_fb(s_,:)] =solve_fb(Para,s_);
end

if and(sSize<3 ,sigma>0)
    
   ssPol=[.7 .7 0 -0.9];
    get_root_ss_nag= @(num,ssPol,user,iflag) getSteadyState(num,ssPol,Para,user,iflag)  ;
    [ssPol,~,exitflag]=c05qb(get_root_ss_nag,ssPol,'xtol',1e-10);
    n_ss=ssPol(1:2);
    phiss=ssPol(3)
    xss=ssPol(4)
end

Para.n_fb=n_fb;
Para.c_fb=c_fb;
Para.x_fb=min(x_fb(1,:));
xMin=min(x_fb(1,:));
xMax=-xMin*.5;
% Define the functional space
for s=1:sSize
V(s) = fundefn(ApproxMethod,OrderOfApprx ,xMin,xMax,orderspli);
end
xGrid=linspace(xMin,xMax,xGridSize)';
xGridSize=length(xGrid);
Para.xMax=xMax;
Para.xMin=xMin;
Para.xGrid=xGrid;
%% Initialize the value function with the deterministic guess
gTrue=Para.g;
g=sum(Para.g.*pi(1,:))*ones(1,sSize);
Para.g=g;
s_=1;





%[n_fb_det(s_,:),c_fb_det(s_,:),x_fb_det(s_,:)] =solve_fb(Para,s_);


for xind=1:xGridSize
    x=xGrid(xind);

for s_=1:sSize
 %   if ~(x>x_fb_det(s_))
%xprime0=ones(1,sSize)*min(x_fb_det(s_));
%n(xind,s_,:)=n_fb_det;
%   else
xprime0=x*ones(1,sSize);
n(xind,s_,:)=(fzero(@(n) res_imp_det(n,x) ,n_fb(1)))*ones(1,sSize);
xprime(xind,s_,:)=xprime0;    
 %   end
    u(xind,s_,:)=util(squeeze(n(xind,s_,:))');
Eu(xind,s_)=sum(pi(s_,:).*squeeze(u(xind,s_,:))');
 
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
Para.g=g;
% Initialize with the no-transfers case
%AMSSNoTransfers=load('~/Golosov-Sargent/Data/temp/AMSS_Solution_no_tr_2shocks.mat');

%for s=1:sSize
%coeff0(s,:)=funfitxy(V(s),xGrid,funeval(AMSSNoTransfers.coeff(s,:)',AMSSNoTransfers.V(s),xGrid ));
%end


%% VALUE FUNCTION ITERATION
options = optimset('Display','off','TolX', 1e-10,'TolCon',1e-10,'TolFun',1e-10);
Para.options=options;
coeff=coeff0;
iter=1;
error=1;
while( error>error_tol && iter<NumIter)
tic
s_=1;

parfor xind=1:xGridSize
    x=xGrid(xind);
    %[ntemp,xprimetemp,c,VNewtemp,exitflagtemp] =solveInnerOpt(x,s_,coeff,V,Para,squeeze(n(xind,s_,:))');
  %if exitflagtemp==0
  %  n(xind,s_,:)=ntemp;
  %  xprime(xind,s_,:)=xprimetemp;
  %  exitflag(xind,s_,:)=exitflagtemp;
  %  VNew(xind,s_)=VNewtemp;
  %end
%    if iter>no_transfer_max_iter
   [n(xind,s_,:),xprime(xind,s_,:),c,VNew(xind,s_),exitflag(xind,s_,:)] =...
        solveInnerOptUsingConsOptFminCon(x,s_,coeff,V,Para,[squeeze(n(xind,s_,:))',squeeze(xprime(xind,s_,:))']);
 %   end
end

for s_=2:sSize
  n(:,s_,:)=n(:,1,:);
  xprime(:,s_,:)=xprime(:,1,:);
  VNew(:,s_)=VNew(:,1);
  exitflag(:,s_,:)=exitflag(:,1,:);
end
%% Update the coeff
coeffold=coeff;
exitflag
for s=1:sSize
%coeffNew(s,:)=funfitxy(V(s),xGrid(logical((exitflag(:,s_,:)==solveflag))),VNew(logical((exitflag(:,s_,:)==solveflag))));
%coeffNew(s,:)=funfitxy(V(s),xGrid,VNew(:,s));
[ coeffNew(s,:)] = FitConcaveValueFunction(V(s),VNew(:,s),xGrid,xGrid(1:2:end)); 
%   if iter>no_transfer_max_iter

%[ coeffNew(s,:)] = SmoothPasting(V(s),VNew(:,s),xGrid,[min(x_fb(1,:))],[0 ],[]);
%    end 
       
ValueDiff(s,:)=funeval(coeffold(s,:)',V(s),xGrid)-funeval(coeffNew(s,:)',V(s),xGrid);
end
coeff=coeffold*(1-grelax)+coeffNew*grelax;
Err(iter,:)=((sum(ValueDiff.^2,2)).^(.5))';
error=Err(iter)
iter=iter+1
toc
save(['~/Golosov-Sargent/Data/temp/' storefilename '.mat'],'coeff','V','Para','n','xprime','Err','xGrid')
end


 
save(['~/Golosov-Sargent/Data/temp/' storefilename '.mat'],'coeff','V','Para','n','xprime','Err','xGrid')




    
