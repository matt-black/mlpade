function [ R32 ] = mitlef_r32 ( alpha, beta, x )
%MITLEF_R32 2nd order approximant R^{3,2}_{\alpha,\beta}
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
