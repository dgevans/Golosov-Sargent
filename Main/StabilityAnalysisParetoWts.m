
 clc
 clear all
 close all
 SetPath

 %{ 
This file solves the G-S economy with BGP preferences of the form
 psi.c^(1-sigma)/1-sigma+(1-psi)log[1-l] with following calibrations

1.] The ratio of productivities is 3.3 and the low productivity is
 normalized to 1
2.] psi is choosen to get a average FE of about 0.5
3.] pareto wts are such that the no-shock problem gives a tax rate of about
    20 percent
4.] Government expenditures are about 11 and 13 percent of the output
5.] beta =0.9
%}

% - XXXXXXX ANMOL - CURRENTLY IT USES A LEGACY METHOD FOR GETTING THE
% BASELINE PARAMETERS. NOW THAT WE HAVE A STEADY STATE CODE WE CAN USE THIS
% TO TARGET SOME AGREED MOMENTS IN OBSERVABLES 
SetParaStruc
theta_1=3; % theta high
theta_2=0;  % theta low
g_l_y=.11; % g low
g_h_y=.13; % g high
n1=1;  
n2=1;
tau=.2;
Y=3
g_Y=mean([g_l_y g_h_y]);

% BASELINE GOVERNMENT EXPENDITURE LEVELS
g=g_Y*Y;


beta=.9;

% BASELINE PARETO WTS
alpha_1=0.69;
alpha_2=1-alpha_1;
Para.n1=n1;
Para.n2=n2;
alpha_1=alpha_1*Para.n1;
alpha_2=alpha_2*Para.n2;

% BASELINE PROBABILITY MATRIX
NewPh=.5;
Para.P=[1-NewPh NewPh;1-NewPh NewPh];

% POPULATE THE PARA STRUC WITH THE BASELINE VALUES
Para.beta=.9;
Para.alpha_1=alpha_1;
Para.alpha_2=alpha_2;
Para.g=[g_l_y g_h_y]*Y;
Para.theta_1=theta_1;
Para.theta_2=theta_2;
Para.btild_1=0;
Para.alpha_1=alpha_1;
Para.alpha_2=alpha_2;
sigma=1
gamma=1
Para.sigma=sigma;
Para.gamma=gamma;


gridSize=5;
alphaMin=.3;
alphaMax=.8;
alphaGrid=linspace(alphaMin,alphaMax,gridSize);
Para.U = @(c,l) UCRRA(c,l,Para);
%ParamGrid=cartprod(alphaGrid,theta1Grid,psiGrid);
ParamGrid=cartprod(alphaGrid);
for i=1:length(ParamGrid)
    tic
    Para.alpha_1=ParamGrid(i,1);
    Para.alpha_2=1-ParamGrid(i,1);
    Para.U = @(c,l) UCRRA(c,l,Para);
    [ x,R,PolicyRule ] = findSteadyState( 0,3,Para);
    %[ A XSS, B, BS ] = LinearApproximation( Para);
    PolicyRules(i,:)=PolicyRule;
    
    C1(i,:)=PolicyRule(1:2);
    C2(i,:)=PolicyRule(3:4);
    L1(i,:)=PolicyRule(5:6);
    L2(i,:)=PolicyRule(7:8);
    
    IC1(i,:)= -(C1(i,:)-L1(i,:).^(1+gamma).*C1(i,:).^(sigma));
    IC2(i,:)= -(C2(i,:)-L2(i,:).^(1+gamma).*C2(i,:).^(sigma));
    
Tau(i,:)=1-L1(i,:).^(1+gamma)./(theta_1*C1(i,:).^(-sigma))    
    Int(i)=C2(i,1)^(-1)./(sum(C2(i,:).^(-1).*Para.P(1,:)));
    X(i)=x;
    RR(i)=R;
   % StabTest(i,:)=max(abs(eigs(BS-eye(4))'));
    toc
end
subplot(1,4,1)
plot(1-ParamGrid,X)
subplot(1,4,2)
plot(1-ParamGrid,RR)
subplot(1,4,3)
plot(1-ParamGrid,PolicyRules(:,5))
subplot(1,4,4)
plot(1-ParamGrid,PolicyRules(:,3))



figure()
plot(1-ParamGrid, Int)


figure()
plot(1-ParamGrid, Tau)


alpha_2=.05

%C=@(rho)  (1+rho)./((1-alpha_2)*(1-rho)+2*rho^2*alpha_2)

%fplot(C,[0,10])