function [SimData]=runSimulation(x0,s0,Para,coeff,V,xGrid,nstore,rhist0,NumSim,flagtr)
g=Para.g;

der_u_c=Para.der_u_c;
der_u_cn=Para.der_u_cn;
der_u_n=Para.der_u_n;
der_u_nn=Para.der_u_nn;


pi=Para.pi;
xHist=zeros(NumSim,1);
sHist=zeros(NumSim,1);
bHist=zeros(NumSim,1);
cHist=zeros(NumSim,1);
nHist=zeros(NumSim,1);

tauHist=zeros(NumSim,1);
% 
% bHist(1)=b_;
% [n0guess,~]=GetInitialApproxPolicy((1/1-g(s0))*b_,xGrid,squeeze(nstore(:,s0,:)));
% 
% get_root_labor0_nag= @(num,n00,user,iflag) getResLaborTime0(num,b_,user,s0,n00,coeff,V,Para,iflag);
% [n00,~,~]=c05qb(get_root_labor0_nag,n0guess,'xtol',1e-10);
% u_c00=der_u_c(n00);
% u_n00=der_u_n(n00);
% 
% xprime00=u_n00.*n00-(n00-g).*u_c00 + u_c00*b_;
% 
%  
%  nHist(1)=n00(s0);
% c00=n00-g(s0);
% xprime0=xprime00(s0);
% sHist(1)=s0;
% xHist(1)=xprime0;
% cHist(1)=c00;
% tau00=1-u_n00./u_c00;
% tauHist(1)=tau00(s0);
tic
sHist(1)=s0;
xHist(1)=x0;


for i_sim=2:NumSim
    x=xHist(i_sim-1);
    s_=sHist(i_sim-1);
[nguess,~]=GetInitialApproxPolicy(x,xGrid,squeeze(nstore(:,s_,:)));  
[n,xprime,c,~,exitflag] =solveInnerOpt(x,s_,coeff,V,Para,nguess);
if flagtr==1

    [n,xprime,c,~,~] =solveInnerOptUsingConsOptFminCon(x,s_,coeff,V,Para,[n,xprime]); 
end

u_c=der_u_c(n);
u_n=der_u_n(n);
Euc=sum(pi(s_,:).*u_c);

   
    tau=1-u_n./u_c;
    
if rhist0(i_sim)<pi(s_,1);
    sHist(i_sim)=1;
else
    sHist(i_sim)=2;
end

xHist(i_sim)=xprime(sHist(i_sim));
nHist(i_sim)=n(sHist(i_sim));
cHist(i_sim)=c(sHist(i_sim));
bHist(i_sim)=xHist(i_sim)*u_c(sHist(i_sim));
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
save('~/Golosov-Sargent/Data/temp/AMSSSimData.mat','SimData')
end
end
        SimData.sHist=sHist;
SimData.xHist=xHist;
SimData.tauHist=tauHist;
SimData.nHist=nHist;
SimData.cHist=cHist;
SimData.bHist=bHist;
save('~/Golosov-Sargent/Data/temp/AMSSSimData.mat','SimData')
