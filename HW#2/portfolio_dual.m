function [obj] = portfolio_dual(mu, Sigma, lambda)

yalmip clear
options = sdpsettings('verbose', 1, 'dualize', 0, 'solver', 'mosek');
rng default

theta = sdpvar(1,1);
beta = sdpvar(length(mu),1);

obj = theta + 0.25*(1/lambda)*((1-lambda)*mu-theta+beta)'*inv(Sigma)*((1-lambda)*mu-theta+beta);

constraints = {};
constraints{end+1} = beta >= 0;

optimize([constraints{:}], obj, options);

theta = double(theta);
beta = double(beta);
obj = theta + 0.25*(1/lambda)*((1-lambda)*mu-theta+beta)'*inv(Sigma)*((1-lambda)*mu-theta+beta);
end