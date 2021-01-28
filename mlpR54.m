function [ R54 ] = mlpR54 ( alpha, beta, x )
% MLPR54 R^{5,4}_{\alpha,\beta}(x)
    pq = pqR54 (alpha, beta);
    R54 = computeRab4o (alpha, beta, pq, x);
end