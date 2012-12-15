function [ Para,V] = BuidGrid( Para)
% This  function defines the grid  and defines the value function
% We have two alternatives. First the user could input either the u2btild
% or Rgrid. this should supersede any other option. Otherwise we use the
% steady state computation and use the default DeltaX,DeltaR paramters to
% set the deviation from the steady state. 

%Para.flagSetu2BtildGrid sets the flag for either using the default grid
%(value =0) or using the user defined grid (Value =1)

%Para.flagSetRdGrid sets the flag for either using the default grid
%(value =0) or using the user defined grid (Value =1) 



% FIND STEADY STATE
[ xSS,RSS,~] = findSteadyState( 0,3,Para);

% CHECK THE FLAG
if isfield(Para,'flagSetu2BtildGrid')
    flagSetu2BtildGrid=Para.flagSetu2BtildGrid;
else
    flagSetu2BtildGrid=0;
end

% USER DEFINED GRID ENDPOINTS
if flagSetu2BtildGrid==1
    disp('Msg :using user defined grid on x')
u2btildMin=Para.u2btildMin;
u2btildMax=Para.u2btildMax;
% DEFAULT GRID ENDPOINTS
else
    disp('Msg :using default grid around SS')

u2btildMin=xSS-Para.DeltaX;
u2btildMax=xSS+Para.DeltaX;

end

% UNIFORMLY SPACE THE GRIDPOINTS
u2btildGrid=linspace(u2btildMin,u2btildMax,Para.u2btildGridSize);
% UPDATE THE PARA STRUCT
Para.u2bdiffGrid=u2btildGrid;
Para.u2btildLL=u2btildMin;
Para.u2btildUL=u2btildMax;



% CHECK THE FLAG

if isfield(Para,'flagSetRGrid')
    flagSetRGrid=Para.flagSetRGrid;
else
    flagSetRGrid=0;
end

% DEFAULT GRID ENDPOINTS
if flagSetRGrid==0
disp('Msg :using default grid around SS')
    RMin=RSS-Para.DeltaR;
    RMax=RSS+Para.DeltaR;

else
    % USER DEFINED GRID ENDPOINTS
    disp('setting RGrid with user inputs')
    RMin=Para.RMin;
    RMax=Para.RMax;
end
% UNIFORMLY SPACE THE GRIDPOINTS
RGrid=linspace(RMin,RMax,Para.RGridSize);
Para.RGrid=RGrid;
GridSize=Para.u2btildGridSize*Para.RGridSize*Para.sSize;

% UPDATE PARATRUC
Para.GridSize=GridSize;
Para.u2btildMin=u2btildMin;
Para.u2btildMax=u2btildMax;
Para.RMax=RMax;
Para.RMin=RMin;

%% DEFINE FUNCTIONAL SPACE
% This uses the CompEcon Library routine `fundefn' to create a functional
% space which stores key settings for the basis polynomials, domain and
% nodes. V(s) is the functional space for the value function given the
% discrete shock

V(1) = fundefn(Para.ApproxMethod,[Para.OrderOfAppx_u2btild Para.OrderOfApprx_R ] ,[u2btildMin RMin],[u2btildMax RMax]);
V(2) = V(1); % 

end

