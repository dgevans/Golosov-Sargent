function [SimData]=runSimulationUsingPolicyRules(b_,s0,Para,coeffN,N,coeff,V,rhist0,NumSim)
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


for i_sim=2:NumSim
    x=xHist(i_sim-1);
    s_=sHist(i_sim-1);
    if rhist0(i_sim)<pi(s_,1);
    sHist(i_sim)=1;
else
    sHist(i_sim)=2;
    end
n(1)=funeval(coeffN(1,:)',N(1),x);
n(2)=funeval(coeffN(2,:)',N(2),x);
c=n-g;
R=1./(beta*(c).*sum(pi(s_,:).*(1./c)));
xprime=(n./(1-n))*(1-psi)+ x.*R-psi;
xHist(i_sim)=xprime(sHist(i_sim));
nHist(i_sim)=n(sHist(i_sim));
if mod(i_sim,1000)==0 || i_sim==NumSim-1
        disp('Running Simulation, t=')
        disp(i_sim)
SimData.sHist=sHist;
SimData.xHist=xHist;
save('~/Golosov-Sargent/Data/temp/AMSSSimDataPol.mat','SimData')
end
end
SimData.sHist=sHist;
SimData.xHist=xHist;
save('~/Golosov-Sargent/Data/temp/AMSSSimDataPol.mat','SimData')
