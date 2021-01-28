function [ R32 ] = mitlef_r32a ( alpha, x )
%MITLEF_R32a 2nd order approximant R^{3,2}_{\alpha,\alpha}
%
    denom = gamma (1+alpha) + ...
        (2 .* gamma(1-alpha).^2) ./ (gamma(1-2.*alpha)) .* x + ...
        gamma (1-alpha) .* x.^2;
    R32 = alpha ./ denom;
end