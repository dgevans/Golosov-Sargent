 
% No transfer - 2 shock case
 MainAMSS(300,'AMSS_ces_no_tr_2_shocks',[.1 .2],0,1,2,.92)
 AMSS_ces_no_tr_2_shocks=load('~/Golosov-Sargent/Data/temp/AMSS_ces_no_tr_2_shocks.mat');
 MainAMSS(300,'AMSS_ces_no_tr_3_shocks',[.1 .15 .2],0,1,2,.92)
 AMSS_ces_no_tr_3_shocks=load('~/Golosov-Sargent/Data/temp/AMSS_ces_no_tr_3_shocks.mat');
 %  transfer - 2 shock case
 MainAMSS(1,'AMSS_ces_tr_2_shocks',[.1 .2],0,1,2,.92)
 AMSS_ces_tr_2_shocks=load('~/Golosov-Sargent/Data/temp/AMSS_ces_tr_2_shocks.mat');
 MainAMSS(1,'AMSS_ces_tr_3_shocks',[.1 .15 .2],0,1,2,.92)
 AMSS_ces_tr_3_shocks=load('~/Golosov-Sargent/Data/temp/AMSS_ces_tr_3_shocks.mat');
 
NumSim=30000;
rhist0=rand(NumSim,1);
s0=1

   
    if(matlabpool('size') == 0)
        matlabpool open 4;
        
    end
    
    
    
    x0=0;



 spmd
    switch labindex
        case 1
            flagtr=0;
[lab1]=runSimulation(x0,s0,AMSS_ces_no_tr_2_shocks.Para,AMSS_ces_no_tr_2_shocks.coeff,AMSS_ces_no_tr_2_shocks.V,AMSS_ces_no_tr_2_shocks.xGrid,AMSS_ces_no_tr_2_shocks.n,rhist0,NumSim,flagtr);
        case 2
            flagtr=0;
[lab2]=runSimulation(x0,s0,AMSS_ces_no_tr_3_shocks.Para,AMSS_ces_no_tr_3_shocks.coeff,AMSS_ces_no_tr_3_shocks.V,AMSS_ces_no_tr_3_shocks.xGrid,AMSS_ces_no_tr_3_shocks.n,rhist0,NumSim,flagtr);
        case 3
            flagtr=1;
[lab3]=runSimulation(x0,s0,AMSS_ces_tr_2_shocks.Para,AMSS_ces_tr_2_shocks.coeff,AMSS_ces_tr_2_shocks.V,AMSS_ces_tr_2_shocks.xGrid,AMSS_ces_tr_2_shocks.n,rhist0,NumSim,flagtr);
        case 4
            flagtr=1;
[lab4]=runSimulation(x0,s0,AMSS_ces_tr_3_shocks.Para,AMSS_ces_tr_3_shocks.coeff,AMSS_ces_tr_3_shocks.V,AMSS_ces_tr_3_shocks.xGrid,AMSS_ces_tr_3_shocks.n,rhist0,NumSim,flagtr);           
    end
end

AMSS_ces_no_tr_2_shocks.SimData=lab1{1};
AMSS_ces_no_tr_3_shocks.SimData=lab2{2};
AMSS_ces_tr_2_shocks.SimData=lab3{3};
AMSS_ces_tr_3_shocks.SimData=lab4{4};

save('~/Golosov-Sargent/Data/temp/AMSS_ces_no_tr_2_shocks.mat','AMSS_ces_no_tr_2_shocks')
save('~/Golosov-Sargent/Data/temp/AMSS_ces_no_tr_3_shocks.mat','AMSS_ces_no_tr_3_shocks')
save('~/Golosov-Sargent/Data/temp/AMSS_ces_tr_2_shocks.mat','AMSS_ces_tr_2_shocks')
save('~/Golosov-Sargent/Data/temp/AMSS_ces_tr_3_shocks.mat','AMSS_ces_tr_3_shocks')

