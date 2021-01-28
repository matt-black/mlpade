function [ R72 ] = mlpR72 ( alpha, beta, x )
% MLPR72 R^{7,2}_{\alpha,\beta}(x)
% 
% only valid for alpha < beta
% if \alpha = \beta use `mitlef_r72a`
    pq = pqR72 (alpha, beta);
    R72 = computeRab4o (alpha, beta, pq, x);
end
