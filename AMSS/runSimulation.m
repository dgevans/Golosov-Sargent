function runSimulation(b_,s0,Para,coeff,V,NumSim)
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
tau0=1-((1-psi)/psi)*(n00-g(s0))/(1-n00);
tauHist(1)=tau0;
nguess=[nDet(xprime00) nDet(xprime00) ];



for i_sim=2:NumSim
    x=xHist(i_sim-1);
    s_=sHist(i_sim-1);
  get_root_labor_nag= @(num,n,user,iflag) getResLaborFsolve(num,n,x,s_,coeff,V,Para,user,iflag)  ;
 
 [n,~,exitflag]=c05qb(get_root_labor_nag,nguess,'xtol',1e-10);
    
%[n,~,exitflag]=fsolve(@(n)getResLaborFsolve(x,s_,n,coeff,V,Para),nguess ,options);
c=n-g;
R=1./(beta*(c).*sum(pi(s_,:).*(1./c)));
% Use Implementability to get xprime
xprime=(n./(1-n))*(1-psi)+ x.*R-psi;
flagConsBind=0;
for s=1:sSize
if xprime(s)>xMax
    xprime(s)=xMax;
    flagConsBind=1;
elseif  xprime(s)<xMin
     xprime(s)=xMin;
    flagConsBind=1;
end
end
if flagConsBind==1
[n,~,exitflag,~,~]=fsolve(@(n) getResInitialLaborFsolve(x,s_,n,Para,xprime),n,options);
c=n-g;
R=1./(beta*(c).*sum(pi(s_,:).*(1./c)));
end
tau=1-((1-psi)/psi).*(n-g)./(1-n);
if rand<pi(s_,1);
    sHist(i_sim)=1;
else
    sHist(i_sim)=2;
end
xHist(i_sim)=xprime(sHist(i_sim));
nHist(i_sim)=n(sHist(i_sim));
bHist(i_sim)=xHist(i_sim)*psi/(c(sHist(i_sim)));
tauHist(i_sim)=tau(sHist(i_sim));
if exitflag==1
nguess=n;
end

end
% Plot Simulations
figure()
subplot(2,2,1)
   X.Data=nHist;
    X.sHist=sHist;
    X.name={'$n$'};
    PlotSimul(X,1);
    title('Labor')
subplot(2,2,2)
   X.Data=bHist;
    X.sHist=sHist;
    X.name={'$b$'};
    PlotSimul(X,1)
    title('Assets');
    subplot(2,2,3)
   X.Data=tauHist;
    X.sHist=sHist;
    X.name={'$\tau$'};
    
    PlotSimul(X,1);
    title('Taxes');
    subplot(2,2,4)
   X.Data=xHist;
    X.sHist=sHist;
    X.name={'$x$'};
    PlotSimul(X,1);
    title('x=u_c b');
    print(gcf,'-dpng',[plotpath 'FigAMSSSimulation.png'])
    disp(tauHist)
   end
