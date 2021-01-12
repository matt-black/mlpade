function [ R72 ] = mitlef_r72 ( alpha, beta, x )
%MITLEF_R72 4th order approximant R^{7,2}_{\alpha,\beta}
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
         -gbma/gbp2a];
    % solve system
    pq = A \ b;                         % [p1,p2,p3,q0,q1,q2,q3]
    % compute
    R72 = (1 / gbma) .* ...
          (pq(1) + pq(2).*x + pq(3).*x.^2 + x.^3) ./ ...
          (pq(4) + pq(5).*x + pq(6).*x.^2 + pq(7).*x.^3 + x.^4);
end
