function [ R63 ] = mlpR63a ( alpha, x )
% MLPR63 R_{\alpha,\alpha}^{6,3}(x)
    pq = pqR63a (alpha);
    R63 = computeRa4o (alpha, pq, x);
end