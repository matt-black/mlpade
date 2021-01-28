function [ R32 ] = mlpR32a ( alpha, x )
% MLPR32A R^{3,2}_{\alpha,\alpha}(x)
%
    denom = gamma (1+alpha) + ...
        (2 .* gamma(1-alpha).^2) ./ (gamma(1-2.*alpha)) .* x + ...
        gamma (1-alpha) .* x.^2;
    R32 = alpha ./ denom;
end