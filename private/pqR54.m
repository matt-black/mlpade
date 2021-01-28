function [ pq ] = pqR54 ( alpha, beta )
% PQR54 solve linear system for R^{5,4}
%
% OUTPUTS: [p1,p2,p3,q0,q1,q2,q3]

    % precompute some terms
    r1d = - gamma (beta-alpha) / gamma (beta);
    r2d = gamma (beta-alpha) / gamma (beta+alpha);
    r3d = - gamma (beta-alpha) / gamma (beta+2*alpha);
    r4d = gamma (beta-alpha) / gamma (beta+3*alpha);
    m2d = gamma (beta-alpha) / gamma (beta-2*alpha);
    m3d = gamma (beta-alpha) / gamma (beta-3*alpha);
    % specify system
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
    % solve
    pq = A \ b;
end