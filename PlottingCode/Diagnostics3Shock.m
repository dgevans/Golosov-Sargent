clear all
clc
close all
rootDirTemp=pwd;
rootDir=rootDirTemp(1:end-length(['/PlottingCode'])); % get root directory
run([rootDir '/Main/SetPath.m'])

ex(1).casename='sigmalow'; 
ex(2).casename='sigmaMed';
ex(3).casename='sigmaHigh'

i=1
% load the coeffecients
BellmanData=load(['Data/temp/c' ex(i).casename '.mat']);
 BellmanData.Para.plotpath=[rootDir sl 'Graphs/3Shocks'];
 BellmanData.Para.datapath=[rootDir sl 'Data/temp/'];

 PlotValueFunction(BellmanData.Para)