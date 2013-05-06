% Run simulations
clear all
clc
close all
AMSSTransfers2Shocks=load('~/Golosov-Sargent/Data/temp/AMSS_Solution_no_tr_2shocks.mat');
AMSSTransfers3Shocks=load('~/Golosov-Sargent/Data/temp/AMSS_Solution_no_tr_3shocks.mat');
NumSim=30000;
rhist0=rand(NumSim,1);
s0=1
b_=0;
err=[];
try
    matlabpool('size')
catch err
end

if isempty(err)
   
    if(matlabpool('size') == 0)
        matlabpool open 2;
        
    end
    
    
end

%[AMSSNoTransfers.SimDataPol]=runSimulationUsingPolicyRules(b_,s0,AMSSNoTransfers.Para,AMSSNoTransfers.coeffN,AMSSNoTransfers.N,AMSSNoTransfers.coeff,AMSSNoTransfers.V,rhist0,NumSim)
spmd
    switch labindex
        case 1  
[lab1]=runSimulation(0,s0,AMSSTransfers2Shocks.Para,AMSSTransfers2Shocks.coeff,AMSSTransfers2Shocks.V,rhist0,NumSim);
        case 2
[lab2]=runSimulation(0,s0,AMSSTransfers3Shocks.Para,AMSSTransfers3Shocks.coeff,AMSSTransfers3Shocks.V,rhist0,NumSim);
    end
end
AMSSTransfers2Shocks.SimData=lab1{1};
AMSSTransfers3Shocks.SimData=lab2{2};
save('~/Golosov-Sargent/Data/temp/AMSSTransfers2Shocks.mat','AMSSTransfers2Shocks');
save('~/Golosov-Sargent/Data/temp/AMSSTransfers3Shocks.mat','AMSSTransfers3Shocks');
