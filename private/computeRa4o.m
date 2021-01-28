function [ R ] = computeRa4o ( alpha, pq, x )
% COMPUTERA4O
    R = (-1/gamma(-alpha)) .* ...
        (pq(1) + pq(2).*x + x.^2) ./ ...
        (pq(3) + pq(4).*x + pq(5).*x.^2 + pq(6).*x.^3 + x.^4);
end

% COMPUTE BY PARTIAL FRACTION
% has problems b/c of weird floating point errors such that you don't 
% actually get roots that are complex conjugates
% function [ R ] = computeRa4o ( pq, x )
%     [c, r, ~] = residue ([1, pq(2), pq(1)], ...
%         [1, pq(6), pq(5), pq(4), pq(3)]);
%     R = c(1)./(x-r(1)) + c(2)./(x-r(2)) + c(3)./(x-r(3)) + c(4)./(x-r(4));
% end
