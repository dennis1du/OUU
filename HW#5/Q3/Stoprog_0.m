function [optimal] = Stoprog(param)

% parameters
p_0 = 4;
p_i = 14;
p_j = 24;
[I,J] = size(param.cost);
S = size(param.outcome,2);

% initialization
timerVal = tic;
options = sdpsettings('verbose', 0, 'solver', 'mosek');

% variables
x = sdpvar(I,J, 'full');
q = sdpvar(I,1);
t = sdpvar(J,1);
y = sdpvar(J,S, 'full');

% objective function
obj = p_i * sum(q) + sum(sum(param.cost.*x)) - p_j * sum(t) + (p_0 + p_j)...
    * sum(sum(param.prob .* y));

% constraints
Cons = [x(:) >= 0, q >= 0, t >= 0, y(:) >= 0];
Cons = [Cons, sum(x,2) == q, sum(x,1) == t', q <= param.cap];
Cons = [Cons, y >= repelem(t,1,S) - param.outcome];

% solve
solvesdp(Cons, obj, options);

% output
optimal.obj = double(obj);
optimal.q = double(q);
optimal.t = toc(timerVal);

end