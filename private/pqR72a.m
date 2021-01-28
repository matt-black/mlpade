function [ pq ] = pqR72a ( alpha )
        % precompute some fractions
    r1d =  gamma (-alpha) / gamma (alpha);
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
         0 1  0   0   0   -1];
    b = [0;
         0;
         -1;
         0;
         gamma(-alpha)/gamma(-3*alpha);
         -gamma(-alpha)/gamma(-2*alpha)];
    % solve system
    pq = A \ b;
end