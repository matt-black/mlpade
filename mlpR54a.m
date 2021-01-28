function [ R54 ] = mlpR54a ( alpha, x )
% MLPR54A R^{5,4}_{\alpha,\alpha}(x)
    
    r1d = gamma (-alpha) / gamma (alpha);
    r2d = -gamma (-alpha) / gamma(2*alpha);
    r3d = gamma (-alpha) / gamma(3*alpha);
    m2d = -gamma(-alpha) / gamma (-2*alpha);
    % formulate linear system
    A = [1 0 r1d  0   0  0;
         0 1 r2d r1d  0  0;
         0 0 r3d r2d r1d 0;
         0 0  0  -1  m2d 0;
         1 0  0   0   -1 -m2d;
         0 1  0   0   0  -1];
    b = [0;
         0;
         -1;
         0;
         0;
         m2d];
    pq = A \ b;
    R54 = computeRa4o (alpha, pq, x);
end