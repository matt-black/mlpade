function [ R63 ] = mlpR63 ( alpha, beta, x )
% MLPR63 R_{\alpha,\beta}^{6,3}(x)
    r1d = - gamma(beta-alpha) / gamma (beta);
    r2d = gamma (beta-alpha) / gamma (beta+alpha);
    r3d = - gamma (beta-alpha) / gamma (beta+2*alpha);
    r4d = gamma (beta-alpha) / gamma (beta+3*alpha);
    r5d = - gamma (beta-alpha) / gamma (beta+4*alpha);
    m1d = gamma (beta-alpha) / gamma (beta-2*alpha);
    % specify A, b
    A = [1 0 0 r1d  0   0   0;
         0 1 0 r2d r1d  0   0;
         0 0 1 r3d r2d r1d  0;
         0 0 0 r4d r3d r2d r1d;
         0 0 0 r5d r4d r3d r2d;
         0 1 0  0   0  -1  m1d;
         0 0 1  0   0   0   -1];
    b = [0;
         0;
         0;
         -1;
         -r1d;
         gamma(beta-alpha)/gamma(beta-3*alpha);
         -m1d];
    % solve the system
    pq = A \ b;
    R63 = computeRab4o (alpha, beta, pq, x);
end