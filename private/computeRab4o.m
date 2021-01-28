function [ R ] = computeRab4o ( alpha, beta, pq, x )
% COMPUTERAB4O compute 4th order pade approximant
    R = (1 / gamma(beta-alpha)) .* ...
        (pq(1) + pq(2).*x + pq(3).*x.^2 + x.^3) ./ ...
        (pq(4) + pq(5).*x + pq(6).*x.^2 + pq(7).*x.^3 + x.^4);
end

% COMPUTE BY PARTIAL FRACTION
% has problems b/c of weird floating point errors such that you don't 
% actually get roots that are complex conjugates
% function [ R ] = computeRab4o ( pq, x )
% % COMPUTERAB4O compute 4th order pade approximant
%     [c, r, k] = residue ([1, pq(3), pq(2), pq(1)], ...
%         [1, pq(7), pq(6), pq(5), pq(4)]);
%     assert (isempty (k), 'mlpade/mlp:computeRab4o :: k wasnt empty');
%     R = c(1)./(x-r(1)) + c(2)./(x-r(2)) + c(3)./(x-r(3)) + c(4)./(x-r(4));
% end