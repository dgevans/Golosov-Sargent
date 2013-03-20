function [ res] = GetCalibrationFrischElasticity (x,AvfFETarget,theta_1,theta_2,tau,g_Y,n1,n2)
gamma=x(1);
Y=x(2);
% GBC with no debt
Trans=((tau-g_Y)/(n1+n2))*Y;
% Labor choices from Golosov's notes
l1Num=(1-tau)*theta_1-gamma*Trans;
l2Num=(1-tau)*theta_2-gamma*Trans;
l1Den=(1+gamma)*(1-tau)*theta_1;
l2Den=(1+gamma)*(1-tau)*theta_2;
l1=l1Num/l1Den;
l2=l2Num/l2Den;
% Average FE u_l/(l u_ll) = 1/l -1 
AvgFE=(n1*(1/l1-1)+n2*(1/l2-1))/(n1+n2);
res(1)=Y-theta_1*l1*n1-theta_2*l2*n2;
res(2)=(AvfFETarget-AvgFE);
end

