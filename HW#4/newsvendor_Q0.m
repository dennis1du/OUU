% initialization
clear all
load('params.mat')
yalmip clear
options = sdpsettings('verbose', 1, 'dualize', 0, 'solver', 'mosek');
rng default

% set variables
N = length(v);
x = sdpvar(N,1);
lambda = sdpvar(1,1);
alpha = sdpvar(1,1);
beta = sdpvar(N,1);
gamma = sdpvar(N,N);
rn = sdpvar(repelem(N,N),repelem(1,N));
sn = sdpvar(repelem(1,N),repelem(1,N));

% objective function
obj = lambda+1/epsilon*(alpha+beta'*mu+trace(gamma'*(Sigma+mu*mu')));

% constraints
Q = zeros(N,N);
sumrn = 0;
sumsn = 0;
for i = 1:N
    sumrn = sumrn + rn{i};
    sumsn = sumsn + sn{i};
end
Con2 = [x >= 0, gamma >= 0];
Con2 = [Con2, [gamma 0.5*beta; 0.5*beta' alpha] >= 0];
Con2 = [Con2, [gamma-Q 0.5*(beta-sumrn); 0.5*(beta-sumrn)' alpha+lambda-(c-v)'*x-sumsn] >= 0];
for i = 1:N
    index = zeros(N,1);
    Con2 = [Con2, [Q 0.5*rn{i}; 0.5*rn{i}' sn{i}] >= 0];
    index(i) = 1;
    Con2 = [Con2, [Q 0.5*(rn{i}+v(i)*index); 0.5*(rn{i}+v(i)*index)' sn{i}-v(i)*x(i)] >= 0];
end

% solve
solvesdp(Con2, obj, options);
obj_approx = double(obj);
