function [ R63 ] = mlpR63a ( alpha, x )
% MLPR63 R_{\alpha,\alpha}^{6,3}(x)
    r1d = gamma(-alpha)/gamma(alpha);
    r2d = -gamma(-alpha)/gamma(2*alpha);
    r3d = gamma(-alpha)/gamma(3*alpha);
    r4d = -gamma(-alpha)/gamma(4*alpha);
    m1d = gamma(-alpha)/gamma(-2*alpha);
    m2d = gamma(-alpha)/gamma(-3*alpha);
    % specify the system
    A = [1 0 r1d  0   0   0;
         0 1 r2d r1d  0   0;
         0 0 r3d r2d r1d  0;
         0 0 r4d r3d r2d r1d;
         1 0  0   0  -1  m1d;
         0 1  0   0   0  -1];
    b = [0;
         0;
         -1;
         0;
         m2d;
         -m1d];
    % solve & compute
    pq = A \ b;
    R63 = computeRa4o (alpha, pq, x);
end