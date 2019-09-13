function [optimal] = Benders(param)

% parameters
p_0 = 4;
p_i = 14;
p_j = 24;
[I,J] = size(param.cost);
S = size(param.outcome,2);
ex_pt = cell(J*S,1);

% initialization
timerVal = tic;
options = sdpsettings('verbose', 0, 'solver', 'mosek');

% variables
x = sdpvar(I,J, 'full');
q = sdpvar(I,1);
t = sdpvar(J,1);
y = sdpvar(J,S, 'full');

% master obj
obj = p_i * sum(q) + sum(sum(param.cost .* x)) - p_j * sum(t) + (p_0 + p_j)...
    * sum(sum(param.prob .* y));

% constraints
Cons = [x(:) >= 0, q >= 0, t >= 0];
Cons = [Cons, sum(x,2) == q, sum(x,1) == t', q <= param.cap];

% algorithm
optimal.obj = [];
while 1
    for j = 1:J
        for s = 1:S
            if ~isempty(ex_pt{(j-1)*S+s})
                Cons = [Cons, y(j,s) >= (t(j) - param.outcome(j,s))...
                    * ex_pt{(j-1)*S+s}(size(ex_pt{(j-1)*S+s},2))];
            end
        end
    end
    
    solvesdp(Cons, obj, options);
    optimal.obj = [optimal.obj, double(obj)];
    
    % extreme points
    flag = 1;
    epsilon = 0.01;
    y_current = double(y);
    for j = 1:J
        for s = 1:S
            t_j = double(t(j));
            if t_j - param.outcome(j,s) <= 0
                pi = 0;
            else
                pi = 1;
            end
            
            % extreme points set
            if (t_j - param.outcome(j,s)) * pi > y_current(j) + epsilon
                flag = 0;
                if isempty(ex_pt{(j-1)*S+s})
                    n = 0;
                else
                    n = size(ex_pt{(j-1)*S+s},2);
                end
                ex_pt{(j-1)*S+s}(n+1) = pi;
            end
        end
    end
    
    if flag == 1
        optimal.q = double(q);
        optimal.t = toc(timerVal);
        break;
    end
end

end