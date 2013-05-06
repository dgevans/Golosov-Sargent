
clear all
figure()
load('/home/anmol/Dropbox/2011RA/FiscalPolicy/OrganizedCode/Golosov-Sargent/Data/temp/BootStrapIneq.mat')
load('/home/anmol/Dropbox/2011RA/FiscalPolicy/OrganizedCode/Golosov-Sargent/Data/temp/cinequality.mat')
t=10
T=length(SD(1).TauHist)
k=110
K=10
T=110
S=2
for j=t:T
    x=SD(1).xHist(j);
    R=SD(1).RHist(j);
    s_=1
    TauCor(j-t+1)=getTaxMoments(x,R,s_,c,V,domain,PolicyRulesStore,Para,S);
AutoCorr(j)=corr(SD(1).TauHist(t:j),SD(i).TauHist(t-1:j-1));    
end
plot(TauCor,AutoCorr)


for i=k:K
    for j=t:T
        
AutoCorr(i,j-t+1)=corr(SD(i).TauHist(j-100:j),SD(i).TauHist(j-101:j-1));
    end

  
end

plot(AutoCorr')
xlabel('Sample Length')
ylabel('Auto corr taxes')
   print(gcf,'-dpng','figAutoCorrIneq.png') 
clear all
K=100
figure()
load('/home/anmol/Dropbox/2011RA/FiscalPolicy/OrganizedCode/Golosov-Sargent/Data/temp/BootStrapProd.mat')

t=10
T=length(SD(1).TauHist)
for i=1:K
    for j=t:T
AutoCorr(i,j-t+1)=corr(SD(i).TauHist(t:j),SD(i).TauHist(t-1:j-1));
    end

end
plot(AutoCorr')
xlabel('Sample Length')
ylabel('Auto corr taxes')
   print(gcf,'-dpng','figAutoCorrProd.png') 
clear all
