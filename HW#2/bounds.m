function [lower, upper] = bounds(t)
if t <= 0
    lower = 1;
    upper = 1;

elseif t >= 10
    lower = 0;
    upper = 0;

else
    yalmip clear
    options = sdpsettings('verbose', 1, 'dualize', 0, 'solver', 'mosek');

    mu1 = 3;
    mu2 = 13;
    s = [1:9];
    c = zeros(9,1);
    
    for i = t:9
        c(i) = 1;
    end
    
    x1 = sdpvar(9,1);
    x2 = sdpvar(9,1);
    
    obj1 = sum(c'*x1);
    obj2 = sum(c'*x2);
    
    con1 = {};
    con1{end+1} = sum(x1) == 1;
    con1{end+1} = s*x1 == mu1;
    con1{end+1} = (s.^2)*x1 == mu2;
    con1{end+1} = x1 >= 0;
    con2 = {};
    con2{end+1} = sum(x2) == 1;
    con2{end+1} = s*x2 == mu1;
    con2{end+1} = (s.^2)*x2 == mu2;
    con2{end+1} = x2 >= 0;

    optimize([con1{:}], obj1, options);
    optimize([con2{:}], -obj2, options);

    lower = double(sum(c'*x1));
    upper = double(sum(c'*x2));
end