function [ Para,V,xhat,Rhat] = BuidGrid( Para)
% This  function defines the grid  and defines the value function
% We have two alternatives. First the user could input either the x
% or Rgrid. this should supersede any other option. Otherwise we use the
% steady state computation and use the default DeltaX,DeltaR paramters to
% set the deviation from the steady state. 

%Para.flagSetxGrid sets the flag for either using the default grid
%(value =0) or using the user defined grid (Value =1)

%Para.flagSetRdGrid sets the flag for either using the default grid
%(value =0) or using the user defined grid (Value =1) 




% CHECK THE FLAG
if isfield(Para,'flagSetxGrid')
    flagSetxGrid=Para.flagSetxGrid;
else
    flagSetxGrid=0;
end

% USER DEFINED GRID ENDPOINTS
if flagSetxGrid==1
    disp('Msg :using user defined grid on x')
    % FIND STEADY STATE

xMin=Para.xMin;
xMax=Para.xMax;
% DEFAULT GRID ENDPOINTS
else
    disp('Msg :using default grid around SS')
[ xSS,RSS,~] = findSteadyState( 0,3,Para);
xMin=xSS-Para.DeltaX;
xMax=xSS+Para.DeltaX;

end

% UNIFORMLY SPACE THE GRIDPOINTS
xGrid=linspace(xMin,xMax,Para.xGridSize);
% UPDATE THE PARA STRUCT
Para.xGrid=xGrid;
Para.xLL=xMin;
Para.xUL=xMax;



% CHECK THE FLAG

if isfield(Para,'flagSetRGrid')
    flagSetRGrid=Para.flagSetRGrid;
else
    flagSetRGrid=0;
end

% DEFAULT GRID ENDPOINTS
if flagSetRGrid==0
disp('Msg :using default grid around SS for x')
% FIND STEADY STATE
[ xSS,RSS,~] = findSteadyState( 0,3,Para);

    RMin=RSS-Para.DeltaR;
    RMax=RSS+Para.DeltaR;

else
    % USER DEFINED GRID ENDPOINTS
    disp('Msg:setting RGrid with user inputs')
    RMin=Para.RMin;
    RMax=Para.RMax;
end
% UNIFORMLY SPACE THE GRIDPOINTS
S = length(Para.P);
RGrid=linspace(RMin,RMax,Para.RGridSize);
Para.RGrid=RGrid;
GridSize=Para.xGridSize*Para.RGridSize*S;

% UPDATE PARATRUC
Para.GridSize=GridSize;
Para.xMin=xMin;
Para.xMax=xMax;
Para.RMax=RMax;
Para.RMin=RMin;

%% DEFINE FUNCTIONAL SPACE
% This uses the CompEcon Library routine `fundefn' to create a functional
% space which stores key settings for the basis polynomials, domain and
% nodes. V(s) is the functional space for the value function given the
% discrete shock
% type of spline
V(1) = fundefn(Para.ApproxMethod,[Para.OrderOfAppx_x Para.OrderOfApprx_R ] ,[xMin RMin],[xMax RMax],Para.order);
xhat(1)=V(1);
Rhat(1)=V(1);

for s=2:S
V(s) = V(1); % 
xhat(s)=V(1);
Rhat(s)=V(1);
end

end

