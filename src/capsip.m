function x = capsip(targetfun, degree, domain)
% Approximates a function through a polynomial Chebyshev approximation
% computed as a semi-infinite programming problem.
% INPUT:
%   targetfun    := function to approximate
%   polynomfun   := approximation function  
%   domain       := cell array that contains the square interval over
% RETURNS: 
%   An array containing the polynomial coefficients.
% The domain's dimension should math the target function's dimension.

pows = pnomialexps(degree);
coeffdim = length(pows);
domdim = length(domain);

% SIP OBJECTIVE FUNCTION
%   f: R^n+1 -> R; where the n+1 element is 't'. 
objfun = @(x) x(coeffdim + 1);

if domdim == 2
    [W1, W2] = meshgrid(domain{1}, domain{2});
    W = {W1, W2};
    pnomial = @pnomial2d;
elseif domdim == 3
    [W1, W2, W3] = meshgrid(domain{1}, domain{2}, domain{3});
    pnomial = @pnomial3d;
end

x0 = zeros(1, coeffdim + 1);        % initial guess, the last element is 't'
ntheta = 2;                         % number of semi-infinite constraints
A = [];                             % equality constraints
b = [];                             % equality constraints 
Aeq = [];                           % inequality constraints
beq = [];                           % inequality constraints
lb = zeros(1, coeffdim + 1) - Inf;  % x = (x1, x2, x3, x4, x5, x6, x7, t) in R^n+1
ub = zeros(1, coeffdim + 1) + Inf;  % x = (x1, x2, x3, x4, x5, x6, x7, t) in R^n+1
[x, fval, exitflag, output, lambda] = fseminf(objfun, x0, ntheta, @seminfcon, A, b, Aeq, beq, lb, ub);
disp(x);

% disp(fval);
% disp(exitflag);
% disp(output);
% disp(lambda);

% SEMI-INFINITE CONSTRAINTS DEFINITION
% Function signature interface is provided by fseminf authors.
% X := (x, t)
% S := an optional sampling interval
function [c, ceq, K1, K2, S] = seminfcon(x, S)
    % non-linear constraints
    c = [];              
    ceq = [];

    % d(w): target function
    d = targetfun(W);

    % F(x, w)
    f = pnomial(x(1:coeffdim), W, pows);
    
    % t
    [cols, rows] = size(W{1});
    t = zeros(cols, rows) + x(coeffdim + 1);

    % g_1(x, w)
    K1 = d - f - t;
    % g_2(x, w)
    K2 = f - d - t;
end
end