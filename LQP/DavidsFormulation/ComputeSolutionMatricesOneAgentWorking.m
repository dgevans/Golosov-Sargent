clc
clear all
close all
LQaproxBig
% Let U=[c1 c2 l1 l2 psi xi]', xi=[g z1 z2]' , lambda=[lambda1 lambda2] 

% FOC : 
% A0 U = A1 mu +A2 xi + A3 lambda- lambda_ /beta
    
F=[
     c1bar              c2bar               -z1bar*l1bar             0; % R
      omegac1           0                   -omegac1*theta1           -c1bar; % foc c1
     0                  omegac2               0                       -c2bar; % foc c2
      -omegac1*theta1    0                   theta1^1*omegac1+omegal1      z1bar*l1bar; % foc l1
  
    ];

A0=[ % W
    0;% R
    sigma*phi1*l1bar-c1bar; %foc c1
    c2bar-sigma*phi2*l2bar %foc c2
    -phi1*l1bar*(1+gamma);% foc l1
    ];


A1=[
    -gbar   z1bar*l1bar                            ;% R
    0       omegac1*eta1                           ; % c1
    0       0                                      ; %c2
    0       omegal1*omegaz1-theta1*eta1*omegac1    ; %l1
    ];


A2=[
    0 0; % R
    1 0; % c1
    0 1;% c2
    0 0;% l1
    ];

% taxes, transfers and output as a function of controls and shocks
P0_tau=(1-taubar)/(taubar)*[-sigma 0 -gamma ];
P1_tau=(1-taubar)/(taubar)*[0 1 ];
P0_T=Tbar^(-1)*[0 c2bar-phi2*l2bar*sigma 0 ];
P1_T=Tbar^(-1)*[0 0];
P0_y=[0 0 z1bar*l1bar]/ybar;
P1_y=[0  z1bar*l1bar]/ybar;
P0=[P0_tau;P0_T;P0_y];
P1=[P1_tau;P1_T;P1_y];
% Primary deficit
G=[-taubar*ybar 2*Tbar -taubar*ybar];

Iu=[1 0 0 0  ;0 1 0 0  ;0 0 1 0 ];


% FOCs Sigma
Sigma0=Iu*inv(F)*A0;
Sigma1=Iu*inv(F)*A1;
Sigma2=Iu*inv(F)*A2;
Ic=[1 0 0 ;0 1 0 ];

% Laws of Motion for Multipliers
SigmaDeltaLambda=inv(Ic*Sigma2)*Ic*Sigma1;
SigmaDeltaMu=(G*P0*Sigma1 +G*P1+gbar*[1 0 ])./(G*P0*Sigma0);

% Final Solution
Sigma0bar=[Sigma0 zeros(size(Sigma0)) zeros(size(Sigma0))]+[zeros(length(Sigma2),1) (1-beta^-1)*Sigma2];
Sigma1bar=(Sigma1-Sigma2*inv(Ic*Sigma2)*Ic*Sigma1);
Sigma2bar=Sigma2*inv(Ic*Sigma2)*Ic*Sigma1;

% Decomposition
DirectEffect=P0*Sigma1bar+P1;
MultiplierEffect=P0*Sigma0*SigmaDeltaMu*(1-beta)+P0*Sigma2*SigmaDeltaLambda*beta;




%simulation

T=100;
gg=.3775/.3194-1;
pr_g=.04;
ineq_g=.03;
P=[.5 .5];
gh=gg;
gl=-gg;
g=[gl;gh];
z1pr=[pr_g;-pr_g];
z1ineq=[ineq_g;-ineq_g];
AA=[1-beta*P(1) -beta*P(2);-beta*P(1) 1-beta*P(2)];


% type of experiment
% gshocks
z1=[0;0];
g=[gl;gh];
% ineq
%z1=z1ineq;
%g=[0;0];
% productivity
%z1=z1pr;
%g=[0;0];
% annuity values
fg=inv(AA)*g*(1-beta);
fz1=inv(AA)*z1*(1-beta);
fgbar=mean(fg);
fz1bar=mean(fz1);
shocks=[g';z1';];
annuity_shocks=[fg';fz1';];
innovations_f=[fg'-fgbar;fz1'-fz1bar;];



mu0=lambda(1)*0;
lambda0=lambda(end-1:end)*0;
mu(1)=mu0;
lambda_(:,1)=lambda0;
s(1)=1;

for n=1:5
rhist=rand(T,1);

    for t=2:T
        if rhist(t)<P(1)
        s(t)=1;
        else
        s(t)=2;
        end
        % law for mu
    mu(t)=mu(t-1)+SigmaDeltaMu*innovations_f(:,s(t));
    % law for lambda
    lambda_(:,t)=lambda_(:,t-1)-SigmaDeltaLambda*(shocks(:,s(t))-annuity_shocks(:,s(t)));    
    % controls from FOC
    U(:,t)=Sigma0*mu(t)+Sigma1*(shocks(:,s(t)))+Sigma2*(lambda_(:,t)-beta^(-1)*lambda_(:,t-1));
    % star problem - no implementability or Euler equations
    Ustar(:,t)=Sigma1*(shocks(:,s(t)));
    % starstar problem - no euler equations
    Ustarstar(:,t)=Sigma0*mu(t)+Sigma1*(shocks(:,s(t)));
    
    
    TaxesTrasnfers(:,t)=P0*U(:,t)+P1*(shocks(:,s(t))); 
    TaxesTrasnfersStar(:,t)=P0*Ustar(:,t)+P1*(shocks(:,s(t)));
    TaxesTrasnfersStarStar(:,t)=P0*Ustarstar(:,t)+P1*(shocks(:,s(t)));
    
    shocksHist(:,t)=(shocks(:,s(t)));
    end
Output=TaxesTrasnfers(3,2:end);
Taxes=TaxesTrasnfers(1,2:end);
TaxRevenues=(TaxesTrasnfers(1,2:end)+TaxesTrasnfers(3,2:end))*taubar*ybar;
Transfers=TaxesTrasnfers(2,2:end);
OutputStar=TaxesTrasnfersStar(3,2:end);
TaxesStar=TaxesTrasnfersStar(1,2:end);
TaxRevenuesStar=(TaxesTrasnfersStar(1,2:end)+TaxesTrasnfersStar(3,2:end))*taubar*ybar;
TransfersStar=TaxesTrasnfersStar(2,2:end);

OutputStarStar=TaxesTrasnfersStarStar(3,2:end);
TaxesStarStar=TaxesTrasnfersStarStar(1,2:end);
TaxRevenuesStarStar=(TaxesTrasnfersStarStar(1,2:end)+TaxesTrasnfersStarStar(3,2:end))*taubar*ybar;
TransfersStarStar=TaxesTrasnfersStarStar(2,2:end);


RelativeChange=Transfers./TaxRevenues;
Consumption(n,:)=(U(1,2:end));
MRN(n)=mean(RelativeChange);
CorrTauY(n)=corr(TaxesTrasnfers(1,2:end)',TaxesTrasnfers(3,2:end)');
std_tau(n)=std(TaxesTrasnfers(1,2:end)');
end


figure()
subplot(3,1,1)
plot(MRN)
subplot(3,1,2)
plot(CorrTauY)
subplot(3,1,3)
plot(std_tau)





figure()
plot(Output,'k','LineWidth',2)
hold on
plot(Transfers,'r','LineWidth',2)
hold on
plot(Taxes,'b','LineWidth',2)

legend('output','transfers','taxes')

 bbT{1}=0:1;
Inx=find(s(1:T)<2);
    for i=2:length(Inx)-1
        bbT{i}=Inx(i)-1:Inx(i);
    end

ShadePlotForEmpahsis( bbT,'k',.05);  

%%print(gcf,'-dpdf','Simulations.pdf')

figure()
plot(OutputStar,'k','LineWidth',2)
hold on
plot(TransfersStar,'r','LineWidth',2)
hold on
plot(TaxesStar,'b','LineWidth',2)

legend('output*','transfers*','taxes*')

 bbT{1}=0:1;
Inx=find(s(1:T)<2);
    for i=2:length(Inx)-1
        bbT{i}=Inx(i)-1:Inx(i);
    end

ShadePlotForEmpahsis( bbT,'r',.05);  

%print(gcf,'-dpng','SimulationsStar.png')


figure()
plot(OutputStarStar,'k','LineWidth',2)
hold on
plot(TransfersStarStar,'r','LineWidth',2)
hold on
plot(TaxesStar,'b','LineWidth',2)

legend('output**','transfers**','taxes**')

 bbT{1}=0:1;
Inx=find(s(1:T)<2);
    for i=2:length(Inx)-1
        bbT{i}=Inx(i)-1:Inx(i);
    end

ShadePlotForEmpahsis( bbT,'r',.05);  

%print(gcf,'-dpng','SimulationsStarStar.png')

figure()
plot(TaxRevenues,'k','LineWidth',2)
hold on
plot(Transfers,'r','LineWidth',2)
legend('TaxRevenues','transfers')

%ShadePlotForEmpahsis( bbT,'r',.05);  
