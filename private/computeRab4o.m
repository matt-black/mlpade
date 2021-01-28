function [ R ] = computeRab4o ( alpha, beta, pq, x )
% COMPUTERAB4O compute 4th order pade approximant
    gbma  = gamma (beta - alpha);
    R = (1 / gamma(beta-alpha)) .* ...
        (pq(1) + pq(2).*x + pq(3).*x.^2 + x.^3) ./ ...
        (pq(4) + pq(5).*x + pq(6).*x.^2 + pq(7).*x.^3 + x.^4);
end