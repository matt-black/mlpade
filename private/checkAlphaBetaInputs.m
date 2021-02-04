function [ alpha, beta ] = checkAlphaBetaInputs ( alpha, beta )
% CHECKALPHABETAINPUTS make sure specified values are valid
    if nargin == 1
        beta = alpha;
    end
    if isempty (beta)
        beta = alpha
    end
    % validate values with assertions
    assert (alpha > 0 & alpha <= 1, 'mlpade/checkAlphaBetaInputs :: alpha must be in range (0,1]');
    assert (alpha <= beta, 'mlpade/imlp :: alpha must be less than beta')
end

