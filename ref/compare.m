addpath (pwd, '..')
VAL = 500; N = 100;
warning ('off', 'all');

[alpha,beta] = meshgrid (linspace (1e-7, 2, N), ...
                         linspace (1e-7, 2, N));

mlfA = nan (size (N, N)); mlfT = nan (size (N, N));
mlA  = nan (size (N, N)); mlT  = nan (size (N, N));
r72A = nan (size (N, N)); r72T = nan (size (N, N));
r32A = nan (size (N, N)); r32T = nan (size (N, N));
for ai = 1:N
    for bi = 1:N
        % only valid for alpha < beta
        if alpha(ai) > beta(ai), continue, end
        tic;
        % reference 1
        mlfA(ai,bi) = mittag_leffler (alpha(ai), beta(bi), VAL);
        mlfT(ai,bi) = toc; tic;
        % reference 2
        mlA(ai,bi)  = ml (VAL, alpha(ai), beta(bi));
        mlT(ai,bi) = toc; tic;
        % pade approximation
        r72A(ai,bi) = mitlef_r72 (alpha(ai), beta(bi), VAL);
        r72T(ai,bi) = toc;
        % pade approximation
        r32A(ai,bi) = mitlef_r32 (alpha(ai), beta(bi), VAL);
        r32T(ai,bi) = toc;
    end
end