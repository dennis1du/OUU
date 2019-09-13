function [optimal] = Stoprog(param)

% initialization
timerVal = tic;
options = sdpsettings('verbose', 0, 'solver', 'mosek');

% parameters
p_product = 14;
p_sell = 24;
p_deposit = 4;
[I,J] = size(param.cost);
S = size(param.outcome, 2);

% variables
x = sdpvar(I,J, 'full');
y = sdpvar(J,S, 'full');

% constraints
Cons = [x(:) >= 0];
Cons = [Cons, sum(x, 2) <= param.cap];
Cons = [Cons, y(:) >= 0, y >= repelem(sum(x, 1)',1, S) - param.outcome];

% objective function
obj = (p_product - p_sell) * sum(x(:)) + sum(sum(param.cost .* x))...
    + (p_deposit + p_sell) * sum(sum(param.prob .* y));

% solve
solvesdp(Cons, obj, options);

% output
optimal.obj = double(obj);
optimal.q = double(sum(x, 2));
optimal.t = toc(timerVal);

end