function e = mittag_leffler(alpha, beta, z, rho)
% MITTAG_LEFFLER   Generalized Mittag-Leffler function.
%   E = mittag_leffler(alpha, beta, z) evaluates the Mittag-Leffler 
%   function for each element of z with parameters alpha and beta, where
%   alpha, beta are real scalars and alpha > 0.
%
%   E = mittag_leffler(alpha, beta, z, rho)  evaluates the Mittag-Leffler 
%   function with accuracy rho. Other parameters are same as before. The 
%   default value of rho is 10^(-12).
%
%   The Generalized Mittag-Leffler function is defined as
%   E_{alpha, beta}(z) = sum k from 0 to inf of Z^k/gamma(alpha*k + beta)
%
%   Reference: 
%       (1) Rudolf Gorenflo, Yuri Luchko, and Joulia Loutchko
%           Computation of the Mittag-Leffler function and its derivatives.
%           Fractional Calculus & Applied Analysis. 2002, 5(4):12-15.
%       (2) MLF (alpha,beta,Z,P) written by Igor Podlubny, Martin Kacenak.
%
%   Author: Yang Zongze
%   E-mail: yangzognze@gmail.com
%   Date: 2016-12-14

    if nargin <3
        error('mittag_leffler::Not enough input arguments.');
    end
    if nargin < 4
        rho = 10^(-12);
    end
    if length(alpha(:)) ~= 1 || length(beta(:)) ~= 1
        error('mittag_leffler::alpha and beta must be scalars.');
    end
    if alpha <= 0 || ~isreal(beta)
        error('mittag_leffler::alpha must be positive and beta must be real.');
    end
    
    size_of_z = size(z);  % save dimensions of z.
    z = z(:);
    if alpha > 1
        % m = floor(alpha) + 1;
        m = floor(alpha/2) + 1 + (alpha < 2);
        new_z = bsxfun(@(x, y)x.^(1/m)*exp(1i*2*pi*y/m), z, 0:m-1);
        e = sum(mittag_leffler(alpha/m, beta, new_z, rho/m), 2)/m; % Need rho/m?
        e = reshape(e, size_of_z);
        return
    end
    
    q = 1; % How to choose q? q < 1 to ensure lower the time cost.
    large = 10+5*alpha; % 10 + 5*alpha; Can be larger.
    tol = eps;     % Used in equal compare.
    integral_tol = rho; % Used in integral.

    index_z_small = abs(z) < q - tol & abs(z) > 0;
    index_z_large = abs(z) > large;
    index_z_one = (abs(z) >= q - tol & abs(z) <= q + tol);
    index_z_special = abs(z) > q + tol & abs(z) <= large;
    z_absz_small = z(index_z_small);
    z_absz_large = z(index_z_large);
    z_absz_one = z(index_z_one);
    z_absz_special = z(index_z_special);

    e = zeros(size(z));
    e(abs(z) == 0) = 1/gamma(beta);
    
    % for z small
    if ~isempty(z_absz_small)
        k_zero = ceil(max(abs((1 - beta)/alpha), max(log(rho*(1-abs(z_absz_small)))./log(abs(z_absz_small)))));
        % make sure gamma do not overflow. Also this can make k_zero small for z close to 1.
        k_zero = min([k_zero, floor((170-beta)/alpha)]); 
        e_absz_small = sum(bsxfun(@(x,y)x.^y./gamma(beta+alpha*y), z_absz_small, 0:k_zero), 2);
        e(index_z_small) = e_absz_small;
    end
    
    % for z large
    if ~isempty(z_absz_large)
        k_zero = floor(-log(rho)./log(min(abs(z_absz_large))));
        e_absz_large = - sum(bsxfun(@(x,y)x.^y./gamma(beta+alpha*y), z_absz_large, -(1:k_zero)), 2);
        index = abs(angle(z_absz_large)) < alpha*pi/4 + min(pi, alpha*pi)/2;
        e_absz_large(index) = e_absz_large(index) + ...
            1/alpha*z_absz_large(index).^((1-beta)/alpha).*exp(z_absz_large(index).^(1/alpha));
        e(index_z_large) = e_absz_large;
    end

    % for z one
    if ~isempty(z_absz_one)
        k_zero = min([floor((170-beta)/alpha),1500]);
        e_absz_one = sum(bsxfun(@(x,y)x.^y./gamma(beta+alpha*y), z_absz_one, 0:k_zero), 2);
        e(index_z_one) = e_absz_one;
    end
    
    % for z special
    if isempty(z_absz_special)
        e = real(e);
        e = reshape(e, size_of_z);
        return
    end
    
    % beta is large
    % k = max{0, (2z-beta)/alpha, (log(rho)+gammaln(beta))/ln(z)
    if beta > 15 && alpha > 1/2
        k_max = min([floor((170-beta)/alpha), log(realmax)/log(max(abs(z_absz_special))),1500]);
        e(index_z_special) = sum(bsxfun(@(x,y)x.^y./gamma(beta+alpha*y), z_absz_special, 0:k_max), 2);
        e = real(e);
        e = reshape(e, size_of_z);
        return
    end
    
    e_absz_special = zeros(size(z_absz_special));
    index_argz_lt_pi_alpha = abs(angle(z_absz_special)) < pi*alpha - tol;
    index_argz_gt_pi_alpha = abs(angle(z_absz_special)) > pi*alpha + tol;
    index_argz_eq_pi_alpha = ~(index_argz_lt_pi_alpha | index_argz_gt_pi_alpha);
    z_argz_gt_pi_alpha = z_absz_special(index_argz_gt_pi_alpha);
    z_argz_lt_pi_alpha = z_absz_special(index_argz_lt_pi_alpha);
    z_argz_eq_pi_alpha = z_absz_special(index_argz_eq_pi_alpha);

    if beta >= 0
        r_zero = max([ 1, 2*max(abs(z_absz_special)), (-log(pi*rho/6))^alpha]);
    else
        r_zero = max(max((-beta+1)^alpha, 2*max(abs(z_absz_special))),...
            (-2*log(pi*rho/(6*(-beta+2)*(-2*beta)^(-beta))))^alpha);
    end
    
    % Choose numerical integral method.
    nflag = 1;  % Use gauss.
    if nflag == 0
        numerical_integral = @numerical_integral_gauss;
    else
    	numerical_integral = @numerical_integral_romberg;
    end
       
    % z_argz_gt_pi_alpha
    if isempty(z_argz_gt_pi_alpha)
        e_argz_gt_pi_alpha = z_argz_gt_pi_alpha;
    elseif beta <= 1
        e_argz_gt_pi_alpha = ...        % (18)
            numerical_integral(@K, 0, r_zero, integral_tol, alpha, beta, z_argz_gt_pi_alpha);
    else
        epsilon = 1; % How to choose epsilon in (17)?
        e_argz_gt_pi_alpha = ...
            numerical_integral(@K, epsilon, r_zero, integral_tol, alpha, beta, z_argz_gt_pi_alpha)...
            + numerical_integral(@P, -pi*alpha, pi*alpha, integral_tol, alpha, beta, epsilon, z_argz_gt_pi_alpha);
    end

    % z_argz_lt_pi_alpha
    if isempty(z_argz_lt_pi_alpha)
        e_argz_lt_pi_alpha = z_argz_lt_pi_alpha;
    elseif beta <=1  % beta < 1 + alpha
        e_argz_lt_pi_alpha = ...        % (26)
            numerical_integral(@K, 0, r_zero, integral_tol, alpha, beta, z_argz_lt_pi_alpha)...
             + 1./alpha*z_argz_lt_pi_alpha.^((1-beta)/alpha).*exp(z_argz_lt_pi_alpha.^(1/alpha));
    else
        epsilon = min(abs(z_argz_lt_pi_alpha))/2; % + 1/2 - 1e-3;  %TODO how to choose this epsilon?
        e_argz_lt_pi_alpha = ...        % (25)
            numerical_integral(@K, epsilon, r_zero, integral_tol, alpha, beta, z_argz_lt_pi_alpha)...
            + numerical_integral(@P, -pi*alpha, pi*alpha, integral_tol, alpha, beta, epsilon, z_argz_lt_pi_alpha)...
            + 1./alpha*z_argz_lt_pi_alpha.^((1-beta)/alpha).*exp(z_argz_lt_pi_alpha.^(1/alpha));
    end

    % z_argz_eq_pi_alpha
    % delta to ensure condition psilon > |z| in (24)
    if ~isempty(z_argz_eq_pi_alpha)
        % epsilon = (max(abs(z_argz_eq_pi_alpha)) + 1)/2; % How to choose epsilon in (24)? 
        epsilon = max(abs(z_argz_eq_pi_alpha)) + 1/2;
        e_argz_eq_pi_alpha = ...
            numerical_integral(@K, epsilon, r_zero, integral_tol, alpha, beta, z_argz_eq_pi_alpha)...
            + numerical_integral(@P, -pi*alpha, pi*alpha, integral_tol, alpha, beta, epsilon, z_argz_eq_pi_alpha);
        e_absz_special(index_argz_eq_pi_alpha) = e_argz_eq_pi_alpha;
    end
    
    % assemble
    e_absz_special(index_argz_gt_pi_alpha) = e_argz_gt_pi_alpha;
    e_absz_special(index_argz_lt_pi_alpha) = e_argz_lt_pi_alpha;

    e(index_z_special) = e_absz_special;
    e = real(e);
    e = reshape(e, size_of_z);
end

function I = numerical_integral_romberg(fun, lower, upper, tol, varargin)
    m = 6;
    h = (upper - lower)/2^m;
    M = zeros(2^m+1, m+1);
    for i = 0:m
        M((0:2^(m-i):2^m) + 1, i+1) = h*2^(m-i);
    end
    M(1, :) = M(1, :)/2;
    M(end, :) = M(end, :)/2;
    
    r = lower + h*(0:2^m);
    f = fun(r, varargin{:});
    T = f * M;
    k = 1;
    A1 = spdiags([-1/(4^k-1)*ones(m-k+1, 1) 4^k/(4^k-1)*ones(m-k+1, 1)], [0, -1], m-k+2, m-k+1);
    S = T*A1;
    k = 2;
    A2 = spdiags([-1/(4^k-1)*ones(m-k+1, 1) 4^k/(4^k-1)*ones(m-k+1, 1)], [0, -1], m-k+2, m-k+1);
    C = S*A2;
    k = 3;
    A3 = spdiags([-1/(4^k-1)*ones(m-k+1, 1) 4^k/(4^k-1)*ones(m-k+1, 1)], [0, -1], m-k+2, m-k+1);
    R = C*A3;
    if max(abs(R(:, end) - R(:, end-1))) < tol
        I = R(:, end);
        return
    end
    
    Tn = T(:, end);
    Sn = S(:, end);
    Cn = C(:, end);
    Rn = R(:, end);
    
    for k = (m+1):15
        h = h/2;
        T2n = Tn/2 + h*sum(fun(lower+(2*(1:2^(k-1)) - 1)*h, varargin{:}), 2);
        S2n = 4*T2n/3 - Tn/3;
        C2n = 16*S2n/15 - Sn/15;
        R2n = 64*C2n/63 - Cn/63;
        if max(abs(R2n - Rn)./abs(R2n)) < tol
            I = R2n;
            return
        else
            Tn = T2n;
            Sn = S2n;
            Cn = C2n;
            oldRn = Rn;
            Rn = R2n;
        end
    end
    I = R2n;
    warning('numerical_integral_romberg reached to max iteration!!! rel error = %g\n', max(abs(R2n-oldRn)./abs(R2n)))
end

function f = K(r, alpha, beta, z)
    f = bsxfun(@(x, y)K_inner(x, alpha, beta, y), r, z);
end

function f = P(phi, alpha, beta, epsilon, z)
    f = bsxfun(@(x, y)P_inner(x, alpha, beta, epsilon, y), phi, z);
end

function f = K_inner(r, alpha, beta, z)
    if length(r) == 1 && r == 0
        f = zeros(size(z));
        return
    end
    f = r.^((1-beta)/alpha).*exp(-r.^(1/alpha)) .* (r*sin(pi*(1-beta))-z*sin(pi*(1-beta+alpha))) ./ ...
        ((alpha*pi).*(r.^2 - 2*r.*z*cos(alpha*pi) + z.^2));
    if length(z) == 1
        f(r == 0) = 0;   % fix compution error when r = 0 and beta = 1.
    end
end

function f = P_inner(phi, alpha, beta, epsilon, z)
    omega = phi*(1+(1-beta)/alpha)+epsilon^(1/alpha)*sin(phi/alpha);
    f = epsilon^(1+(1-beta)/alpha)*exp(epsilon^(1/alpha)*cos(phi/alpha)).*...
        (cos(omega) + 1i*sin(omega))./(2*alpha*pi*(epsilon*exp(1i*phi) - z));
end

function I = numerical_integral_gauss(fun, lower, upper, tol, varargin)
    N = ceil(-log10(tol)*10); % get from tol?
    [x, omega] = jacobi_gauss(0, 0, N);
    x = x(:)';
    r = (x+1)/2*(upper-lower) + lower;
    f = fun(r, varargin{:});
    I = f*omega(:)*(upper-lower)/2;
end

function [ x, omega ] = jacobi_gauss( alpha, beta, N )
    % code from https://github.com/lrtfm/Jacobi/blob/master/JacobiGauss.m
    j = 0:N;
    T = 2*j + alpha + beta;
    D0 = (beta^2 - alpha^2)./(T .* (T + 2));
    j = 1:N;
    T = 2*j + alpha + beta;
    D1 = 4 * j .* (j + alpha) .* (j + beta) .* ( j + alpha + beta) ./...
        ((T - 1) .* (T.^2) .*(T+1));
    D1 = D1.^(1/2);
    A = diag(D0) + diag(D1, 1) + diag(D1, -1);
    if alpha == 0 && beta == 0
        A(1, 1) = 0;
    end
    [V, D] = eig(A);
    r = 2^(alpha + beta + 1)*gamma(alpha + 1)* gamma(beta + 1)/gamma(alpha + beta + 2);
    x = diag(D, 0);
    omega = (V(1, :).^2 * r)';
end
