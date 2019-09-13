clear all

% parameters
load tele.mat
[I,J] = size(AA);
K = size(DD,1);
S = 100;
c = [1; 2.9; 1; 2.9; 1; 1; 2.9];
r = zeros(I,S);
d = zeros(K,S);
ws = zeros(S,1);

% simulation
rng default
for j = 1:S
    for i = 1:K
        drand = rand;
        if ismember(i, [1, 5, 6, 7])
            if drand <= 0.05
                d(i,j) = 0;
            elseif (0.05 < drand) && (drand <= 0.25)
                d(i,j) = 1;
            elseif (0.25 < drand) && (drand <= 0.75)
                d(i,j) = 2;
            elseif (0.75 < drand) && (drand <= 0.95)
                d(i,j) = 3;
            else
                d(i,j) = 5;
            end
        else
            if drand <= 0.1
                d(i,j) = 0;
            elseif (0.1 < drand) && (drand <= 0.5)
                d(i,j) = 1;
            elseif (0.5 < drand) && (drand <= 0.9)
                d(i,j) = 2;
            else
                d(i,j) = 3;
            end
        end
    end
end

for j = 1:S
    for i = 1:I
        rrand = rand;
        if ismember(i, [2, 4, 7])
            r(i,j) = 1;
        else
            if rrand <= 0.7
                r(i,j) = 1;
            else
                r(i,j) = 0;
            end
        end
    end
end

% initialization
yalmip clear
options = sdpsettings('verbose', 1, 'dualize', 0, 'solver', 'mosek');

% variables
x = sdpvar(I,S, 'full');
y = sdpvar(J,S, 'full');
b = sdpvar(K,S, 'full');

% WS
for j = 1:S
    obj = sum(b(:,j));

    Cons = [b(:) >= 0, x(:) >= 0, y(:) >= 0];
    Cons = [Cons, c' * x(:,j) <= 20];
    Cons = [Cons, AA * y(:,j) <= r(:,j) .* x(:,j)];
    Cons = [Cons, b(:,j) >= d(:,j) - DD * y(:,j)];
    
    solvesdp(Cons, obj, options);
    
    ws(j) = double(obj);
end

ws_mean = mean(ws);
[h,p,ci,stats] = ttest(ws);