function [ R72 ] = mlpR72a ( alpha, x )
% MLPR72A R^{7,2}_{\alpha,\alpha}(x)
    pq = pqR72a (alpha);
    R72 = computeRa4o (alpha, pq, x);
end
