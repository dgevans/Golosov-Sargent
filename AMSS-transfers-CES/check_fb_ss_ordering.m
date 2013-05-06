clear all
%% Set the technology and preference parameters.
sigma=1; % 
gamma=2;
beta=.92; % time discount factor
g=[.1 .2]; % The vector g is the value of the expenditure shock in each state s
sSize=length(g); % This is the dimension of the markov shock
pi=ones(sSize,sSize)/sSize;
%pi=repmat([.7 .1 .2],sSize,1);
A=eye(sSize);
for s=1:sSize
    for sprime=1:sSize
        if ~(s==sprime)
        A(s,sprime)=-beta*pi(s,sprime);
        else
                   A(s,sprime)=1-beta*pi(s,sprime);
        end 
    end
end

% Populate the utility function and the required derivatives
if sigma >1
util=@(n) (n-g).^(1-sigma)./(1-sigma)-n.^(1+gamma)./(1+gamma); % u(n-g,n)
else
util=@(n) log(n-g)-n.^(1+gamma)./(1+gamma); % u(n-g,n)
end
der_u_c=@(n) (n-g).^(-sigma); % d u /d c (as a function of n)
der_u_n=@(n) n.^(gamma); % d u / d n
der_u_cn=@(n) -sigma*(n-g).^(-sigma-1); % d u /dc^2
der_u_nn=@(n) gamma*n.^(gamma-1); % d u / dn^2
res_imp_det =@(n,x) (n-sum(g.*pi(1,:)))^(-sigma)*(n-sum(g.*pi(1,:)))+x-n^(1+gamma)-x/beta;
% 

%%BGP
%util=@(n) psi*log(n-g)+(1-psi)*log(1-n); % u(n-g,n)
%der_u_c=@(n) psi*(n-g).^(-1); % d u /d c (as a function of n)
%der_u_n=@(n) (1-psi)*(1-n).^(-1); % d u / d n
%der_u_cn=@(n) -sigma*(n-g).^(-2); % d u /dc^2
%der_u_nn=@(n) (1-psi)./(1-n).^2; % d u / dn^2
%res_imp_det =@(n,x) psi*(n-sum(g.*pi(1,:)))^(-1)*(n-sum(g.*pi(1,:)))+x- (1-psi)*(n/(1-n)) -x/beta;


Para.der_u_c=der_u_c;
Para.der_u_n=der_u_n;
Para.util=util;
Para.der_u_cn=der_u_cn;
Para.der_u_nn=der_u_nn;
%% Set other parameters
%Para.psi=psi;
Para.sigma=sigma;
Para.gamma=gamma;
Para.g=g;
Para.pi=pi;
Para.beta=beta;
Para.sSize=sSize;
Para.invA=inv(A);
Para.der_u_c=der_u_c;
Para.der_u_n=der_u_n;
Para.util=util;
Para.u_cn=der_u_cn;
Para.u_nn=der_u_nn;

% We now compute the grid for x. 
for s_=1:sSize
[n_fb(s_,:),c_fb(s_,:),x_fb(s_,:)] =solve_fb(Para,s_);
end
x_fb

   ssPol=[.7 .7 0 -0.9];
    get_root_ss_nag= @(num,ssPol,user,iflag) getSteadyState(num,ssPol,Para,user,iflag)  ;
    [ssPol,~,exitflag]=c05qb(get_root_ss_nag,ssPol,'xtol',1e-10);
    n_ss=ssPol(1:2);
    phi_ss=ssPol(3)
    x_ss=ssPol(4)
