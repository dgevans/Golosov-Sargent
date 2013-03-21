function [ res ] = FOCResidualsN(  X,x,R,Vx,VR,Para )
%FOCRESIDUALS Summary of this function goes here
%   Detailed explanation goes here

P = Para.P(1,:);
S = length(P);
N = Para.N;
x = repmat(x',S,1);
R = repmat(R',S,1);
beta = Para.beta;
alpha = Para.alpha(:)';%ensure row vector
alpha_i = repmat(Para.alpha(1:N-1),S,1);
alpha_N = alpha(N);
if numel(Para.theta) == N
    theta = Para.theta(:)';%ensure row vector
    theta_i = repmat(theta(1:N-1),S,1);
    theta_N = repmat(theta(N),S,N-1);
else
    theta_i = theta(:,1:N-1);
    theta_N = repmat(theta(:,N),1,N-1);
end
n = Para.n(:)';
n_i = repmat(n(1:N-1),S,1);
n_N = n(N);

g = Para.g(:);
U = Para.U;

ci = reshape(X(1:(N-1)*S),S,N-1);
X(1:(N-1)*S) = [];
cN = repmat(reshape(X(1:S),S,1),1,N-1);
X(1:S) = [];
li = reshape(X(1:(N-1)*S),S,N-1);
X(1:(N-1)*S) = [];
lN = repmat(reshape(X(1:S),S,1),1,N-1);
X(1:S) = [];
xprime = reshape(X(1:(N-1)*S),S,N-1);
X(1:S*(N-1)) = [];
Rprime = reshape(X(1:(N-1)*S),S,N-1);
X(1:S*(N-1)) = [];
mu = reshape(X(1:(N-1)*S),S,N-1);
X(1:S*(N-1)) = [];
lambda = repmat(reshape(X(1:N-1),1,N-1),S,1);
X(1:N-1) = [];
phi = reshape(X(1:(N-1)*S),S,N-1);
X(1:(N-1)*S) = [];
xi = repmat(reshape(X(1:S),S,1),1,N-1);
X(1:S) = [];
rho = reshape(X(1:(N-1)*S),S,N-1);
X(1:(N-1)*S) = [];

[~,uci,uli,ucci,ulli] = U(ci,li,Para);
[~,ucN,ulN,uccN,ullN] = U(cN,lN,Para);

Euci = repmat(P * uci,S,1);
EucN = repmat(P * ucN,S,1);
Emu_ucN = repmat(P * (ucN.*mu),S,1);

res = [];


resTemp = (x.*ucN)./(beta.*EucN) - ucN.*(cN-ci)-xprime + Rprime.*li.*uli...
    - lN.*ulN;
res = [res;reshape(resTemp,(N-1)*S,1)];

resTemp = P * ( uci.*( Rprime-R ) );
res = [res; resTemp'];

resTemp = uli.*theta_N.*Rprime - theta_i.*ulN;
res = [res; reshape(resTemp,(N-1)*S,1)];

resTemp = sum(n_i.*(theta_i.*li-ci),2) + n_N.*theta_N(:,1).*(lN(:,1)-cN(:,1)) - g;
res = [res;resTemp];


resTemp = Rprime.*uci-ucN;
res = [res;reshape(resTemp,(N-1)*S,1)];

resTemp = alpha_i.*uci + mu.*ucN + lambda.*ucci.*(Rprime-R) - xi.*n_i + rho.*Rprime.*ucci;
res = [res;reshape(resTemp,(N-1)*S,1)];

resTemp = alpha_N(:,1).*ucN(:,1) + sum( mu.*( x.*uccN./(beta.*EucN) - uccN...
    .*(cN-ci) -ucN) - x.*uccN.*Emu_ucN./(beta*EucN.^2),2) - n_N(:,1).*xi(:,1)...
    -sum( rho.*uccN,2);
res = [res;resTemp];

resTemp = alpha_i.*uli + mu.*Rprime.*( uli+li.*ulli )...
    + phi.*theta_N.*Rprime.*ulli + xi.*theta_i.*n_i;
res = [res;reshape(resTemp,(N-1)*S,1)];

resTemp = alpha_N(:,1).*ulN(:,1) - sum( mu.*( ulN + lN.*ullN ) + phi.*theta_i.*ullN,2)...
    +xi(:,1).*theta_N(:,1).*n_N(:,1);
res = [res;resTemp];

resTemp = beta*Vx - mu;
res = [res;reshape(resTemp,(N-1)*S,1)];

resTemp = beta*VR + mu.*li.*uli + lambda.*uci + phi.*uli.*theta_N + rho.*uci;
res = [res;reshape(resTemp,(N-1)*S,1)];

end

