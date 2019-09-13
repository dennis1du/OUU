function [w,b]=ccp(epsilon)

% data loading
load data.mat;
y = cell2mat(y);
[p, n] = size(y);

% modeling
yalmip clear
options = sdpsettings('verbose', 1, 'dualize', 0, 'solver', 'mosek');

w = sdpvar(p,1);
b = sdpvar(1,1);
flag = binvar(n,1);

obj = 1;

constraints = {};
constraints{end+1} = w <= 100;
constraints{end+1} = w >= -100;
constraints{end+1} = b <= 100;
constraints{end+1} = b >= -100;
constraints{end+1} = (w'*y-b)'.*z >= 1-999*flag;
constraints{end+1} = flag'*ones(n,1)/n <= epsilon;

optimize([constraints{:}], obj, options);

w = double(w);
b = double(b);

for i = 1:n
    if z(i) == 1
        scatter(y(1,i), y(2, i), 'o', 'r');
        hold on
    else
        scatter(y(1,i), y(2, i), 'x', 'b');
        hold on
    end
end
x1 = -0.5:0.01:4;
x2 = (b-w(1)*x1)/w(2);
plot(x1,x2);

end