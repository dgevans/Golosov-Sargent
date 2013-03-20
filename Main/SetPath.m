% SETPATH
% This script makes the code sub-directories available on the current path
% relative to the /Main
global  rootDir sl

% THIS MAKES THE PATH COMPATIBLE ACROSS OPERATING PLATFORM
if strcmp(computer,'PCWIN64')
    sl='\';
    coresize=4;
else
    sl='/';
    coresize=8;
end
% 2. GET ROOT
rootDirTemp=pwd;
rootDir=rootDirTemp(1:end-length([sl 'Main'])); % get root directory


% SET paths for external libraries
compeconpath=[rootDir sl 'compecon2011' sl];
knitropath=[rootDir sl 'knitro' sl];


% UPDATE CODE SUBDIRECTORIES
addpath(genpath(rootDir))
addpath(genpath(compeconpath))
addpath(genpath(knitropath))
rmpath(genpath([rootDir sl '.git']))
rmpath(genpath([rootDir sl 'Auxiliary']))

% PATHS for SAVING, PLOTTING etx
texpath= [rootDir sl 'Tex' sl] ;
plotpath= [rootDir sl 'Graphs' sl] ;
datapath=[rootDir sl 'Data' sl] ;
mkdir(texpath)
mkdir(plotpath)
mkdir(datapath)
