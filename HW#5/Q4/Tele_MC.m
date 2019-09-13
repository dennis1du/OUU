clear all

% parameters
load tele.mat
[I,J] = size(AA);
K = size(DD,1);
S = 1000;
c = [1; 2.9; 1; 2.9; 1; 1; 2.9];
r = zeros(I,S);
d = zeros(K,S);

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
x = sdpvar(I,1);
y = sdpvar(J,S, 'full');
b = sdpvar(K,S, 'full');

% objective function
obj =  sum(sum(b))/S;

% constraints
Cons = [b(:) >= 0, x >= 0, y(:) >= 0];
for j = 1:S
    Cons = [Cons, c' * x <= 20];
    Cons = [Cons, AA * y(:,j) <= r(:,j) .* x];
    Cons = [Cons, b(:,j) >= d(:,j) - DD * y(:,j)];
end

% solve
solvesdp(Cons, obj, options);

MC = double(obj);