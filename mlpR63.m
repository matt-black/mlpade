function [ R63 ] = mlpR63 ( alpha, beta, x )
% MLPR63 R_{\alpha,\beta}^{6,3}(x)
    pq  = pqR63 (alpha, beta);
    R63 = computeRab4o (alpha, beta, pq, x);
end