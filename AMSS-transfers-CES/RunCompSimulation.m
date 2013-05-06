% Run simulations
clear all
clc
close all
AMSSTransfers=load('~/Golosov-Sargent/Data/temp/AMSS_Solution_tr.mat');
AMSSNoTransfers=load('~/Golosov-Sargent/Data/temp/AMSS_Solution.mat');
AMSSNoTransfers.xGrid=AMSSNoTransfers.Para.xGrid;
AMSSNoTransfers.Para.flagTransfers=0;
AMSSTransfers.Para.flagTransfers=1;
NumSim=50000;
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
[lab1]=runSimulation(1.2,s0,AMSSNoTransfers.Para,AMSSNoTransfers.coeff,AMSSNoTransfers.V,rhist0,NumSim);
        case 2
[lab2]=runSimulationWithTransfers(1.2,s0,AMSSTransfers.Para,AMSSTransfers.coeff,AMSSTransfers.V,rhist0,NumSim,AMSSNoTransfers);
       case 3
           
[lab3]=runSimulation(-1,s0,AMSSNoTransfers.Para,AMSSNoTransfers.coeff,AMSSNoTransfers.V,rhist0,NumSim);
        case 4
                        
[lab4]=runSimulationWithTransfers(-1,s0,AMSSTransfers.Para,AMSSTransfers.coeff,AMSSTransfers.V,rhist0,NumSim,AMSSNoTransfers);
    end           
end
AMSSNoTransfers.SimData.bzero=lab1{1};
AMSSTransfers.SimData.bzero=lab2{2};
AMSSNoTransfers.SimData.balt=lab3{3};
AMSSTransfers.SimData.balt=lab4{4};
save('~/Golosov-Sargent/Data/temp/AMSSTransfers.mat','AMSSTransfers');
save('~/Golosov-Sargent/Data/temp/AMSSNoTransfers.mat','AMSSNoTransfers');