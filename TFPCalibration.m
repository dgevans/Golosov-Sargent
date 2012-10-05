%load Data/TFPAnnual.mat
 TFPLevel(1)=1;
for i=2:length(TFPAnnual)
    TFPLevel(i)=TFPLevel(i-1)*(1+TFPAnnual(i)/100);
end
[~,detrendedTFP]=hpfilter(TFPLevel,100)
PositiveTFP=detrendedTFP(find(detrendedTFP>0));
NegativeTFP=detrendedTFP(find(detrendedTFP<0));
exp(mean(PositiveTFP))/exp((mean(NegativeTFP)))
x=fminunc (@(x) ErrorInApprox( x,detrendedTFP),[0 0]);
Ratio=abs(exp(max(x))/exp(min(x)))