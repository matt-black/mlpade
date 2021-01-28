function [ R32 ] = mlpR32 ( alpha, beta, x )
% MLPR32 R^{3,2}_{\alpha,\beta}(x)
%
    % precompute gamma terms
    gbeta = gamma (beta);
    gbpa  = gamma (beta + alpha);
    gbma  = gamma (beta - alpha);
    gbm2a = gamma (beta - 2.*alpha);
    % compute leading term
    c_ab = 1 ./ (gbpa .* gbma - gbeta.^2);
    % compute polynomial terms
    p1 = c_ab .* (gbeta.*gbpa - (gbpa .* gbma.^2) ./ gbm2a);
    q0 = c_ab .* (gbeta.^2 .* gbpa ./ gbma - ...
                  (gbeta .* gbpa .* gbma) ./ gbm2a);
    q1 = c_ab .* (gbeta .* gbpa - gbeta.^2 .* gbma ./ gbm2a);
    % now do the computation
    R32 = (1 ./ gbma) .* ((p1 + x) ./ (q0 + q1.*x + x.^2));
end

% NOTE: this doesn't seem to work in MATLAB due to numerics issues(?)
% compute partial fraction weights/poles
% r1 = (-q1 + sqrt (q1^2 - 4*q0)) / 2;
% r2 = (-q1 - sqrt (q1^2 - 4*q0)) / 2;
% c1 = (p1 - r1) / (r2 - r1);
% c2 = (p1 - r2) / (r1 - r2);
% R32 = 2 .* real (c1 ./ (x - r1));
