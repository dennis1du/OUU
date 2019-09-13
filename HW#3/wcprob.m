function [prob] = wcprob()

yalmip clear
options = sdpsettings('verbose', 1, 'dualize', 0, 'solver', 'mosek');

% set variables
alpha = sdpvar(1,1);
beta = sdpvar(2,1);
gamma = sdpvar(2,1);
theta0 = sdpvar([2 2], [1 1]);
theta1 = sdpvar([1 1 1 1]);
theta2 = sdpvar([2 2 2 2], [1 1 1 1]);
theta3 = sdpvar([2 2 2 2], [1 1 1 1]);

% set parameters
t = 1;
sigma = 0.1;
d{1} = [1; 1];
d{2} = [1; -1];
d{3} = [-1; -1];
d{4} = [-1; 1];

% objective function
obj = alpha-sigma*sum(gamma);

% constraints
Cons = gamma >= 0;
Cons = [Cons, theta0{1} >= 0, theta0{2} >= 0];
Cons = [Cons, alpha <= t];
Cons = [Cons, beta-theta0{1}+theta0{2} == 0];
Cons = [Cons, gamma-theta0{1}-theta0{2} == 0];
for i = 1:4
    Cons = [Cons, theta1{i} >= 0, theta2{i} >= 0, theta3{i} >= 0];
    Cons = [Cons, alpha - theta1{i} <= 0];
    Cons = [Cons, beta+theta1{i}*d{i}-theta2{i}+theta3{i} == 0];
    Cons = [Cons, gamma-theta2{i}-theta3{i} == 0];
end

% solve
solvesdp(Cons, -obj, options);

prob = double(obj);

end