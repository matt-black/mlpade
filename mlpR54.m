function [ R54 ] = mlpR54 ( alpha, beta, x )
% MLPR54 R^{5,4}_{\alpha,\beta}(x)
    r1d = - gamma (beta-alpha) / gamma (beta);
    r2d = gamma (beta-alpha) / gamma (beta+alpha);
    r3d = - gamma (beta-alpha) / gamma (beta+2*alpha);
    r4d = gamma (beta-alpha) / gamma (beta+3*alpha);
    m2d = gamma (beta-alpha) / gamma (beta-2*alpha);
    m3d = gamma (beta-alpha) / gamma (beta-3*alpha);
    
    A = [1 0 0 r1d 0 0 0;
         0 1 0 r2d r1d 0 0;
         0 0 1 r3d r2d r1d 0;
         0 0 0 r4d r3d r2d r1d;
         1 0 0  0   -1 m2d m3d;
         0 1 0  0   0   -1 m2d;
         0 0 1  0   0   0   -1];
    b = [0;
         0;
         0;
         -1;
         -gamma(beta-alpha)/gamma(beta-4*alpha);
         m3d
         -m2d];
    pq = A \ b;
    % compute
    R54 = computeRab4o (alpha, beta, pq, x);
end