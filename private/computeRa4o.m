function [ R ] = computeRa4o ( alpha, pq, x )
% COMPUTERA4O
    R = (-1/gamma(-alpha)) .* ...
        (pq(1) + pq(2).*x + x.^2) ./ ...
        (pq(3) + pq(4).*x + pq(5).*x.^2 + pq(6).*x.^3 + x.^4);
end

