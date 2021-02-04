function [ pq ] = solvePQcoeffs ( alpha, beta, mn )
% SOLVEPQCOEFFS solve the linear system for a given '(m,n)' approximant
%
% INPUTS:
%   alpha : 1st argument to ML function
%   beta  : 2nd arg. to ML function
%   mn    : type of approximant
% OUTPUTS:
%   pq    : vector of polynomial coefficients [p, q]

    if isempty (beta), beta = alpha; end
    % validate inputs
    [alpha, beta] = checkAlphaBetaInputs (alpha, beta);
    switch mn
      case { '54', '(5,4)' }
        if alpha == beta, pq = pqR54a (alpha);
        else, pq = pqR54 (alpha, beta);
        end
      case { '63', '(6,3)' }
        if alpha == beta, pq = pqR63a (alpha);
        else, pq = pqR63 (alpha, beta);
        end
      case { '72', '(7,2)' }
        if alpha == beta, pq = pqR72a (alpha);
        else, pq = pqR72 (alpha, beta);
        end
      otherwise
        error ('mlpade/solvePQcoeffs :: invalid m/n combination')
    end
end