clear all

% parameters
load tele.mat
[I,J] = size(AA);
K = size(DD,1);
S = 100;
c = [1; 2.9; 1; 2.9; 1; 1; 2.9];
r = zeros(I,S);
d = zeros(K,S);
eev = zeros(S,1);

% initialization
yalmip clear
options = sdpsettings('verbose', 1, 'dualize', 0, 'solver', 'mosek');

% solve x_EV
prob1 = [0.05; 0.2; 0.5; 0.2; 0.05];
d1 = [0; 1; 2; 3; 5];
prob2 = [0.1; 0.4; 0.4; 0.1];
d2 = [0; 1; 2; 3];
r_ev = [0.7; 1; 0.7; 1; 0.7; 0.7; 1];
d_ev1 = d1'*prob1;
d_ev2 = d2'*prob2;
d_ev = [d_ev1; d_ev2; d_ev2; d_ev2; d_ev1; d_ev1; d_ev1; d_ev2; d_ev2; d_ev2];

x_ev = sdpvar(I,1);
y_ev = sdpvar(J,1);
b_ev = sdpvar(K,1);

obj_ev =  sum(b_ev);

Cons_ev = [b_ev >= 0, x_ev >= 0, y_ev >= 0];
Cons_ev = [Cons_ev, c' * x_ev <= 20];
Cons_ev = [Cons_ev, AA * y_ev <= r_ev .* x_ev];
Cons_ev = [Cons_ev, b_ev >= d_ev - DD * y_ev];

solvesdp(Cons_ev, obj_ev, options);

x_EV = double(x_ev);

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
            if rand <= 0.7
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
y = sdpvar(J,S, 'full');
b = sdpvar(K,S, 'full');

% WS
for j = 1:S
    obj = sum(b(:,j));

    Cons = [b(:) >= 0, y(:) >= 0];
    Cons = [Cons, AA * y(:,j) <= r(:,j) .* x_EV];
    Cons = [Cons, b(:,j) >= d(:,j) - DD * y(:,j)];
    
    solvesdp(Cons, obj, options);
    
    eev(j) = double(obj);
end

eev_mean = mean(eev);
[h,p,ci,stats] = ttest(eev);