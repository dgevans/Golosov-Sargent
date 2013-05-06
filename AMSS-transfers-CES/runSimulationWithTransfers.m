function [SimData]=runSimulationWithTransfers(b_,s0,Para,coeff,V,rhist0,NumSim)
g=Para.g;
psi=Para.psi;
pi=Para.pi;
beta=Para.beta;
sSize=Para.sSize;
options=Para.options;
plotpath=Para.plotpath;
xMax=Para.xMax;
xMin=Para.xMin;
nDet=@(x) (psi+x*((beta-1)/beta))/(1+x*((beta-1)/beta));
xHist=zeros(NumSim,1);
sHist=zeros(NumSim,1);
bHist=zeros(NumSim,1);
cHist=zeros(NumSim,1);
nHist=zeros(NumSim,1);

tauHist=zeros(NumSim,1);

exitflagsim=zeros(NumSim,1);
bHist(1)=b_;
n0guess=nDet((1/1-g(s0))*b_);
get_root_labor0_nag= @(num,n00,user,iflag) getResLaborTime0(num,b_,user,s0,n00,coeff,V,Para,iflag);
 [n00,~,exitflag0]=c05qb(get_root_labor0_nag,n0guess,'xtol',1e-10);
nHist(1)=n00;
c00=n00-g(s0);
xprime00=(n00/(1-n00))*(1-psi)+ (psi/c00)*b_-psi;
sHist(1)=s0;
xHist(1)=xprime00;
cHist(1)=c00;
tau0=1-((1-psi)/psi)*(n00-g(s0))/(1-n00);
tauHist(1)=tau0;
nguess=[nDet(xprime00) nDet(xprime00) ];
tic

P=Para.pi;
S=length(P(1,:));
CumP=[];

for s_ = 1:S
    CumP(s_,1)=0;
for s=2:S
    CumP(s_,s)=CumP(s_,s-1)+P(s_,s);
end
CumP(s_,S+1)=1;
end


for i_sim=2:NumSim
    x=xHist(i_sim-1);
    s_=sHist(i_sim-1);
  %[nguess,xprimeguess,~,~] =solveInnerOpt(x,s_,AMSSNoTransfers.coeff,AMSSNoTransfers.V,Para);
  [nguess,xprimeguess,~,~] =solveInnerOpt(x,s_,coeff,V,Para);
  [n,xprime,c,~,~] =solveInnerOptUsingConsOptFminCon(x,s_,coeff,V,Para,[nguess,xprimeguess]);
  
sHist(i_sim)=sum(~(CumP(sHist(i_sim-1),:)-rhist0(i_sim)>0));

tau=1-((1-psi)/psi).*(n-g)./(1-n);
xHist(i_sim)=xprime(sHist(i_sim));
nHist(i_sim)=n(sHist(i_sim));
cHist(i_sim)=c(sHist(i_sim));
bHist(i_sim)=xHist(i_sim)*psi/(c(sHist(i_sim)));
tauHist(i_sim)=tau(sHist(i_sim));


    
     if mod(i_sim,1000)==0 || i_sim==NumSim-1
        disp('Running Simulation, t=')
        disp(i_sim)
        toc
        tic
        SimData.sHist=sHist;
SimData.xHist=xHist;
SimData.tauHist=tauHist;
SimData.nHist=nHist;
SimData.cHist=cHist;
SimData.bHist=bHist;
save('~/Golosov-Sargent/Data/temp/AMSSSimData_tr.mat','SimData')
end
end
        SimData.sHist=sHist;
SimData.xHist=xHist;
SimData.tauHist=tauHist;
SimData.nHist=nHist;
SimData.cHist=cHist;
SimData.bHist=bHist;
save('~/Golosov-Sargent/Data/temp/AMSSSimData_tr.mat','SimData')
