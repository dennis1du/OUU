function [optimal] = Benders(param)

% set up
timerVal = tic;
options = sdpsettings('verbose', 0, 'solver', 'mosek');

% parameters
p_product = 14;
p_sell = 24;
p_deposit = 4;
[I,J] = size(param.cost);
S = size(param.outcome, 2);
%delta = 1e-2;
epsilon = 0.01;

% variables
x = sdpvar(I,J, 'full');
gamma = sdpvar(J,S, 'full');
%pi = sdpvar(J,S, 'full');

% Step 0
z_upper = 1e5;
ext_pt = cell(J,S);
%ext_ray = cell(J*S,1);

% initialization
cons_master = [x(:) >= 0];
cons_master = [cons_master, sum(x, 2) <= param.cap];
obj_master = (p_product - p_sell) * sum(x(:)) + sum(sum(param.cost .* x))...
    + (p_deposit + p_sell) * sum(sum(param.prob .* gamma));
optimal.obj = [];

% algorithm
while 1
    % Step 1
    % add extreme point constraints
    for j = 1:J
        for s = 1:S
            if ~isempty(ext_pt{j,s})
               cons_master = [cons_master, (-param.outcome(j,s) +...
            sum(x(:,j))) * ext_pt{j,s}(size(ext_pt{j,s}, 2)) <= gamma(j,s)];
            end
        end
    end    
    
    % solve master problem
    solvesdp(cons_master, obj_master, options);
    optimal.obj = [optimal.obj, double(obj_master)];
    
    % 1st return
    x_hat = double(x);
    gamma_hat = double(gamma);
    z_lower = double(obj_master);
    
    % Step 2
    % solve dual objective function
    for j = 1:J
        for s = 1:S
%             cons_dual = [pi(j,s) >= 0, pi(j,s) <= 1];
%             obj_dual = (-param.outcome(j,s) + sum(x_hat(:,j))) * pi(j,s);
%             solvesdp(cons_dual, -obj_dual, options);
%             z_d(j,s) = double(obj_dual);
%             pi_hat(j,s) = double(pi(j,s));
            if sum(x_hat(:,j)) - param.outcome(j,s) <= 0
                pi = 0;
            else
                pi = 1;
            end
            z_d(j,s) = (sum(x_hat(:,j)) - param.outcome(j,s)) * pi;
            
            % Step 4
            % update ext_pt set
            if z_d(j,s) > gamma_hat(j,s)
                if isempty(ext_pt{j,s})
                    n = 0;
                else
                    n = size(ext_pt{j,s}, 2);
                end
                ext_pt{j,s}(n+1) = pi;
            end
        end
    end
    
    % 2nd return
    z_hat = (p_product - p_sell) * sum(x_hat(:)) + sum(sum(param.cost .* x_hat))...
    + (p_deposit + p_sell) * sum(sum(param.prob .* z_d));
    
    % update
    if z_hat < z_upper
        z_upper = z_hat;
        %x_opt = x_hat;
    end
    
    % Step 3
    % evaluation & output
%    if (z_upper - z_lower <= delta * min(abs(z_upper), abs(z_lower)))
    if z_upper - z_lower <= epsilon
        %optimal.x = x_opt;
        optimal.q = sum(x_hat, 2);
        optimal.t = toc(timerVal);
        break;
    end    
end

end