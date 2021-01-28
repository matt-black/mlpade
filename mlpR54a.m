function [ R54 ] = mlpR54a ( alpha, x )
% MLPR54A R^{5,4}_{\alpha,\alpha}(x)
    pq  = pqR54a (alpha);
    R54 = computeRa4o (alpha, pq, x);
end