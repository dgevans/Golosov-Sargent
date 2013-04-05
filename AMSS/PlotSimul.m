function PlotSimul(X,flagFig)
if nargin==1
figure()
end
T=length(X.Data);

hold on
 bbT{1}=0:1;
Inx=find(X.sHist(1:T)<2);
    for i=2:length(Inx)-1
        bbT{i}=Inx(i)-1:Inx(i);
    end
plot((1:T),X.Data(1:T))
axis tight
ShadePlotForEmpahsis( bbT,'r',.05);  
ylabel(X.name,'Interpreter','Latex')

end
