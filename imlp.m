function [ x ] = imlp ( alpha, beta, y, mn )
% IMLP Pade approximant of inverse Mittag-Leffler function

% INPUTS:
%    alpha : must be in range (0,1]
%    beta  : must be >= alpha
%    x     : value to compute E_{\alpha,\beta} at
%    mn    : type of approximation to compute (R^{m,n})
%            implemented (m,n) pairs are: (3,2), (7,2)
% OUTPUTS:
%    L   :  Pade approximant of E_{\alpha,\beta}(-x)

    if nargin < 4, mn = '(7,2)'; end
    % make sure x inputs are valid (positive)
    if any (y < 0)
        warning ('mlpade/imlp :: negative inputs were automatically flipped to positive values')
        y = abs (y);
    end
    % checks on alpha, beta inputs
    assert (alpha > 0 & alpha <= 1, 'mlpade/imlp :: alpha must be in range (0,1]');
    assert (alpha <= beta, 'mlpade/imlp :: alpha must be less than beta')
    if (alpha == 1 && beta == 1)
        warning ('mlpade/imlp :: alpha=1,beta=1 is just exp')
        x = - log(y);
        return;
    end

    if strcmp (mn, '32') || strcmp (mn, '(3,2)') % 2nd order approximation
        if alpha == beta
            % eqn. (38)
            x = -gamma(1-alpha)./(y.*gamma(1-2*alpha)) + ...
                sqrt (gamma(1-alpha)^2./(gamma(1-2*alpha).^2 .* y.^2) - ...
                      ((1+alpha)/(1-alpha)).*(1-(1./(gamma(alpha).*y))));
        else
            error ('mlpade/imlp :: not implemented')
        end
    else                                % 4th order approximation
        a_eq_b = (alpha == beta);
        % compute `c` coefficient
        if a_eq_b, c = gamma (-alpha);
        else,      c = gamma (beta-alpha);
        end
        % compute p_i, q_i by solving linear system
        switch mn
          case {'54', '(5,4)'}
            if a_eq_b, pq = pqR54a (alpha);
            else,      pq = pqR54 (alpha, beta);
            end
          case {'63', '(6,3)'}
            if a_eq_b, pq = pqR63a (alpha);
            else,      pq = pqR63 (alpha, beta);
            end
          case {'72', '(7,2)'}
            if a_eq_b, pq = pqR72a (alpha);
            else,      pq = pqR72 (alpha, beta);
            end            
          otherwise
            error ('mlpade/imlp :: invalid approximation type, %s', mn);
        end
        % find unique positive root of polynomial for system
        % NOTE: need arrayfun here because `polyRoot` only works for scalar `y` inputs
        x = arrayfun (@(yy) polyRoot (alpha, beta, pq, yy), y);
    end
end


function [ val ] = polyRoot ( alpha, beta, pq, y )
% POLYROOT find unique positive root of polynomial system
% SEE: eqns. 40, 41 of paper
    if alpha == beta
        % eqn. (41)
        c_a = gamma (-alpha);
        p = [c_a*y, c_a*y*pq(end), c_a*y*pq(end-1)+1, ...
             c_a*y*pq(end-2)+pq(2), c_a*y*pq(end-3)+pq(1)];
    else
        % eqn. (40)
        c_ab = gamma (beta-alpha);
        p = [c_ab*y, c_ab*y*pq(end)-1, c_ab*y*pq(end-2)-pq(3), ...
             c_ab*y*pq(end-3)-pq(2), c_ab*y*pq(end-4)-pq(1)];
    end
    % find roots of specified polynomial
    rtz = real (roots (p));
    % inverse is unique positive root
    ind = find (rtz == max (rtz), 1);
    val = rtz(ind);
    % checks that there is a unique positive root
    if sum (sign(rtz) > 0) > 1
        warning ('mlpade/imlp:polyRoot :: system gave multiple positive roots');
    end
    if rtz(ind) < 0
        val = NaN; return
    end
    assert (rtz(ind) > 0, 'mlpade/imlp:polyRoot :: largest root was negative');
end

