% This program solves recrusive AMSS - optimal policy problem for the
% separable BGP case : u(c,1-n)=psi log (c) + (1-psi) log (1-n)
clc
clear all
close all
%% Set the technology and preference parameters.
psi=.2; % psi is the relative weight on consmumption in the utility function
beta=.9; % time discount factor
g=[.1 .3]; % The vector g is the value of the expenditure shock in each state s
sSize=length(g); % This is the dimension of the markov shock
pi=ones(sSize,sSize)/sSize;
%pi=repmat([.7 .1 .2],sSize,1);


deltaG=linspace(0,.4,15)
for i=1:15
g(1)=0;
g(2)=.1+deltaG(i)
n_fb=g+psi*(1-g);
c_fb=psi*(1-g);
uc_fb=1./(1-g);
Euc_fb=sum(pi(1,:).*uc_fb);
int_fb=uc_fb./(beta*Euc_fb);
x_fb(i,:)=(((1-psi)).*(n_fb./(1-n_fb))-psi).*(1-int_fb).^(-1);
end


plot(deltaG,x_fb)