%% PAPER FIGURES
% recreate (roughly) figures from the paper
addpath (pwd, '..')
warning ('off', 'all');

%% Figure 1
% see page 8
abps = {[0.25,1], [0.5,1];
    [0.75,1], [0.95,1]};

xs = linspace (0.01, 10, 25);

figure(1), clf
for sp = 1:4
    alpha = abps{sp}(1);
    beta  = abps{sp}(2);
    ref   = mittag_leffler (alpha, beta, -xs);
    app7   = mlpR32 (alpha, beta, xs);
    subplot (2, 2, sp)
    plot (xs, ref, 'b', xs, app7, 'g^')
    legend ('ml', '(3,2)')
end

%% Figure 2
% see page 8
alphas = [0.25, 0.5];
xs = linspace (0.01, 10, 25);
figure(2), clf
for sp = 1:2
    ref = mittag_leffler (alphas(sp), alphas(sp), -xs);
    app7 = mlpR32a (alphas(sp), xs);
    subplot (1, 2, sp)
    plot (xs, ref, 'b', xs, app7, 'g^')
    if sp == 1
        ylim ([0 0.3])
    else
        ylim ([0 0.55])
    end
    legend ('ml', '(3,2)')
end

%% Figure 3
% see page 15
abps = {[0.25,1], [0.5,1];
    [0.75,1], [0.95,1]};

xs = linspace (.001, 10, 25);

figure(3), clf
for sp = 1:4
    alpha = abps{sp}(1);
    beta  = abps{sp}(2);
    ref   = mittag_leffler (alpha, beta, -xs);
    app5  = mlpR54 (alpha, beta, xs);
    app7  = mlpR72 (alpha, beta, xs);
    subplot (2, 2, sp)
    plot (xs, app5, 'm*', xs, app7, 'ro', xs, ref, 'b')
    legend ('5,4)','(7,2)','ml')
end

%% Figure 4
% see page 16, beta > 1

abps = {[1, 1.01], [1, 2]};

xs = linspace (.001, 10, 25);

figure(4), clf
for sp = 1:2
    alpha = abps{sp}(1);
    beta  = abps{sp}(2);
    ref   = mittag_leffler (alpha, beta, -xs);
    app5  = mlpR54 (alpha, beta, xs);
    app7  = mlpR72 (alpha, beta, xs);
    subplot (1, 2, sp)
    plot (xs, app5, 'm*', xs, app7, 'ro', xs, ref, 'b')
    legend ('(5,4)','(7,2)','ml')
end

%% Figure 5
% see page 16, beta > 1

alphas = [0.5 0.75];

xs = linspace (.001, 10, 25);

figure(5), clf
for sp = 1:2
    ref   = mittag_leffler (alphas(sp), alphas(sp), -xs);
    app5  = mlpR54a (alphas(sp), xs);
    app7  = mlpR72a (alphas(sp), xs);
    subplot (1, 2, sp)
    plot (xs, app5, 'm*', xs, app7, 'ro', xs, ref, 'b')
    legend ('(5,4)','(7,2)','ml')
end

%% Figure 6
abps = {[0.75,1], [0.75,0.75];
    [0.95,1], [1,1.01]};
xs = logspace (-4, 1, 100);

figure(6), clf
for sp = 1:4
    alpha = abps{sp}(1);
    beta  = abps{sp}(2);
    ref   = mittag_leffler (alpha, beta, -xs);
    if alpha == beta
        app7 = mlpR72a (alpha, xs);
        app5 = mlpR54a (alpha, xs);
    else
        app5 = mlpR54 (alpha, beta, xs); 
        app7 = mlpR72 (alpha, beta, xs);
    end
    err5 = abs (ref - app5);
    err7 = abs (ref - app7);
    subplot (2, 2, sp)
    plot (xs, err5, 'b', xs, err7, 'g')
    set (gca, 'YScale', 'log')
    legend ('(5,4)','(7,2)')
    switch sp
        case 1, ylim ([1E-15 1E-3])
        case 2, ylim ([1E-15 1E-1])
        case 3, ylim ([1E-7  1E-1])
        otherwise, ylim ([1E-10 1E-1])
    end
    title (sprintf ('\\alpha=%.2f, \\beta=%.2f', alpha, beta))
end

%% Figure 7
abps = {[0.75,1], [0.5,0.5];
    [0.95,1], [1,2]};
xs = logspace (-4, 1, 100);

figure(7), clf
for sp = 1:4
    alpha = abps{sp}(1);
    beta  = abps{sp}(2);
    ref   = mittag_leffler (alpha, beta, -xs);
    app7 = mlp (alpha, beta, xs, '(7,2)');
    app3 = mlp (alpha, beta, xs, '(3,2)');
    subplot (2, 2, sp)
    plot (xs, app3, 'g^', xs, app7, 'ro', xs, ref, 'b')
    legend ('(3,2)', '(7,2)', 'ml')
    title (sprintf ('\\alpha=%.2f, \\beta=%.2f', alpha, beta))
end

%% Figure 8
abps = {[0.75,1], [0.5,0.5];
    [0.95,1], [1,2]};
xs = logspace (-4, 1, 100);

figure(8), clf
for sp = 1:4
    alpha = abps{sp}(1);
    beta  = abps{sp}(2);
    ref   = mittag_leffler (alpha, beta, -xs);
    if alpha == beta
        app7 = mlpR72a (alpha, xs);
        app3 = mlpR32a (alpha, xs);
    else
        app7 = mlpR72 (alpha, beta, xs);
        app3 = mlpR32 (alpha, beta, xs);
    end
    subplot (2, 2, sp)
    plot (xs, abs (app7-ref), 'b', xs, abs(app3-ref), 'r')
    legend ('(7,2)', '(3,2)')
    title (sprintf ('\\alpha=%.2f, \\beta=%.2f', alpha, beta))
    set (gca, 'YScale', 'log')
end