clear all

% parameters
load tele.mat
[I,J] = size(AA);
K = size(DD,1);
c = [1; 2.9; 1; 2.9; 1; 1; 2.9];
prob1 = [0.05; 0.2; 0.5; 0.2; 0.05];
d1 = [0; 1; 2; 3; 5];
prob2 = [0.1; 0.4; 0.4; 0.1];
d2 = [0; 1; 2; 3];
r = [0.7; 1; 0.7; 1; 0.7; 0.7; 1];
d_ev1 = d1'*prob1;
d_ev2 = d2'*prob2;
d = [d_ev1; d_ev2; d_ev2; d_ev2; d_ev1; d_ev1; d_ev1; d_ev2; d_ev2; d_ev2];

% initialization
yalmip clear
options = sdpsettings('verbose', 1, 'dualize', 0, 'solver', 'mosek');

% variables
x = sdpvar(I,1);
y = sdpvar(J,1);
b = sdpvar(K,1);

% objective function
obj =  sum(b);

% constraints
Cons = [b >= 0, x >= 0, y >= 0];
Cons = [Cons, c' * x <= 20];
Cons = [Cons, AA * y <= r .* x];
Cons = [Cons, b >= d - DD * y];

% solve
solvesdp(Cons, obj, options);

EV = double(obj);
x_EV = double(x);
