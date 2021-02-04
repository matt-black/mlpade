function [ R ] = compute4OrdPoly ( pq, x, varargin )
% COMPUTE4ORDPOLY direct computation of the Pade approximant from a solved system (p,q)
%
% INPUTS:
%   pq              : polynomial coefficient vector (probably from `solvePQcoeffs`)
%   x               : value(s) to compute approximant at
%   alpha           : 1st argument to ML function
%   beta (optional) : 2nd arg. to ML function
% OUTPUTS:
%   R               : value of approximant at `x`
    if (numel (varargin) == 1)
        R = computeRa4o (varargin{1}, pq, x);
    else
        R = computeRab4o (varargin{1}, varargin{2}, pq, x);
    end
end

