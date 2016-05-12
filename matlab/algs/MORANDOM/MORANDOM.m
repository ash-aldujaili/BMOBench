function [pf,ps,fc]= MORANDOM(problem,l , u, numEvaluations, k)
%% FUNCTION MORANDOM
%   A multi-objective optimizer based on uniform random sampling
% problem : function handler problem : R^n -> R^k, gives a vector 1xk
% l : lower bound of the dimension space, 1xn
% u : upper bound of the dimension space, 1xn
% numEvaluations : evaluation budget

%-- Initialize the variables --------------------------------------%
% problem specific 
delta    = u - l;
n        = length(l);
func     = @(x) problem( delta.*x + l);
% get the samples:
c = rand(numEvaluations,n); % centre point of the node, representative state
fc = zeros(numEvaluations,k);
% evaluate them
for j = 1 : numEvaluations
	fc(j,:) = func(c(j,:));
end
% get the approximation set
front = paretofront(fc);
pf = fc(front,:);
ps = bsxfun(@plus, bsxfun(@times,delta, c(front,:)),l);
end
