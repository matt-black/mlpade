function [ R72 ] = mitlef_r72a ( alpha, x )
%MITLEF_R72a 4th order approximant R^{7,2}_{\alpha,\beta} when \alpha = \beta
%
    % precompute some fractions
    r1d =  gamma (-alpha);
    r2d = -gamma (-alpha) / gamma (2*alpha);
    r3d =  gamma (-alpha) / gamma (3*alpha);
    r4d = -gamma (-alpha) / gamma (4*alpha);
    r5d =  gamma (-alpha) / gamma (5*alpha);
    % specify A matrix
    A = [1 0 r1d  0   0   0 ;
         0 1 r2d r1d  0   0 ;
         0 0 r3d r2d r1d  0 ;
         0 0 r4d r3d r2d r1d;
         0 0 r5d r4d r3d r2d;
         0 0  0   0   0   -1];
    b = [0;
         0;
         -1;
         0;
         gamma(-alpha)/gamma(-3*alpha);
         -gamma(-alpha)/gamma(-2*alpha)];
    % solve system
    pq = A \ b;                         % [p2,p3,q0,q1,q2,q3]
    % compute
    R72 = (-1/gamma(-alpha)) .* ...
          (pq(1) + pq(2).*x + x.^2) ./ ...
          (pq(3) + pq(4).*x + pq(5).*x.^2 + pq(6).*x^3 + x.^4);
end
