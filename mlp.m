function [ y ] = mlp ( alpha, beta, x, mn )
% MLP Pade approximant R^{m,n}_{\alpha,\beta}(x) \approx E_{\alpha,\beta}(-x)
%
% INPUTS:
%    alpha : must be in range (0,1]
%    beta  : must be >= alpha
%    x     : value to compute E_{\alpha,\beta} at
%    mn    : type of approximation to compute (R^{m,n})
%            implemented (m,n) pairs are: (3,2), (7,2)
% OUTPUTS:
%    y   :  Pade approximant of E_{\alpha,\beta}(-x)
% 
% unlike all other functions, this one does input sanitation/error checking
% note that this really just dispatches to the underlying functions for the actual computation
    if nargin < 4, mn = '72', end
    % make sure x inputs are valid (positive)
    if any (x < 0)
        warning ('mlpade/mlp :: negative inputs were automatically flipped to positive values')
        x = abs (x);
    end
    % checks on alpha, beta inputs
    assert (alpha > 0 & alpha <= 1, 'mlpade/mlp :: alpha must be in range (0,1]');
    assert (alpha <= beta, 'mlpade/mlp :: alpha must be less than beta')
    if (alpha == 1 && beta == 1)
        warning ('mlpade/mlp :: alpha=1,beta=1 is just exp')
        y = exp (-x);
        return;
    end
    % dispatch to output function
    switch mn
      case {'32', '(3,2)'}
        if alpha == beta
            y = mlpR32a (alpha, x);
        else
            y = mlpR32  (alpha, beta, x);
        end
      case {'72', '(7,2)'}
        if alpha == beta
            y = mlpR72a (alpha, x);
        else
            y = mlpR72  (alpha, beta, x);
        end
      otherwise
        error ('mlpade/mlp :: invalid approximation type, %s', mn);
    end
end