function [ Loss] = ErrorInApprox(x,detrendedTFP)
theta_l=x(1);
theta_h=x(2);
PositiveTFP=detrendedTFP(find(detrendedTFP>0));
NegativeTFP=detrendedTFP(find(detrendedTFP<0));
Loss= ((mean(NegativeTFP)-theta_l)^2*length(NegativeTFP)+ (mean(PositiveTFP)-theta_h)^2*length(PositiveTFP))^0.5;

end

