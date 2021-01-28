function [ R72 ] = mlpR72 ( alpha, beta, x )
% MLPR72 R^{7,2}_{\alpha,\beta}(x)
% 
% only valid for alpha < beta
% if \alpha = \beta use `mitlef_r72a`
    % precompute gamma terms
    gbeta = gamma (beta);
    gbpa  = gamma (beta + alpha);
    gbma  = gamma (beta - alpha);
    gbm2a = gamma (beta - 2*alpha);
    gbp2a = gamma (beta + 2*alpha);
    gbp3a = gamma (beta + 3*alpha);
    gbp4a = gamma (beta + 4*alpha);
    gbp5a = gamma (beta + 5*alpha);
    % precompute some fractions
    r1d = - gbma / gbeta;
    r2d = gbma / gbpa;
    r3d = - gbma / gbp2a;
    r4d = gbma / gbp3a;
    r5d = - gbma / gbp4a;
    r6d = gbma / gbp5a;
    % specify A matrix
    A = [1 0 0 r1d  0   0   0 ;
         0 1 0 r2d r1d  0   0 ;
         0 0 1 r3d r2d r1d  0 ;
         0 0 0 r4d r3d r2d r1d;
         0 0 0 r5d r4d r3d r2d;
         0 0 0 r6d r5d r4d r3d;
         0 0 1  0   0   0   -1];
    b = [0;
         0;
         0;
         -1;
         gbma/gbeta;
         -gbma/gbpa;
         -gbma/gbm2a];
    % solve system
    pq = A \ b;                         % [p1,p2,p3,q0,q1,q2,q3]
    R72 = computeRab4o (alpha, beta, pq, x);
end
