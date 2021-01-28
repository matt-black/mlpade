% figure out how much time is spent on input sanitization for the `mlp`
% function versus the no-checking functions

clear, clc, close all
addpath (pwd, '..')
addpath (pwd, '..', 'fun');

abps = {[0.75,1], [0.5,0.5];
    [0.95,1], [1,2]};
xs = logspace (-4, 1, 5000);

tml7 = nan (4, 2);
tml3 = nan (4, 2);
for sp = 1:4
    alpha = abps{sp}(1);
    beta  = abps{sp}(2);
    % compute from `mlp` function
    tic
    mlp (alpha, beta, xs, '(7,2)');
    tml7(sp,1) = toc; tic
    mlp (alpha, beta, xs, '(3,2)');
    tml3(sp,1) = toc;
    % compute from implementation
    tic
    if alpha == beta
        mlpR72a (alpha, xs);
    else
        mlpR72 (alpha, beta, xs);
    end
    tml7(sp,2) = toc;
    tic
    if alpha == beta
        mlpR32a (alpha, xs);
    else
        mlpR32 (alpha, beta, xs);
    end
    tml3(sp,2) = toc;
end

tml3 = tml3 .* 10^6;
tml7 = tml7 .* 10^6;

fprintf ('AVG overhead R32 : %.2f us\n', mean(tml3(:,1)-tml3(:,2)));
fprintf ('AVG overhead R72 : %.2f us\n', mean(tml7(:,1)-tml7(:,2)));
fprintf ('R72 vs. R32 : %.2f, %.2f us\n', ...
    mean (tml7(:,1)-tml3(:,1)), mean (tml7(:,2)-tml3(:,2)));
