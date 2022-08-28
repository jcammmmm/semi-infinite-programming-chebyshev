d = @(W1, W2) log(W1 + W2).*sin(W1);
capsipf(d, 3);

function x = capsipf(targetfun, degree)
% Approximates a function through a polynomial Chebyshev approximation
% computed as a semi-infinite programming problem.
% INPUT:
%   targetfun    := function to approximate
%   polynomfun   := approximation function  
%   domain       := array that contains the square interval over
% RETURNS: 
%   An array containing the polynomial coefficients.
% The domain's dimension should math the target function's dimension.

pows = getpows(degree);
domdim = length(pows);

% SIP OBJECTIVE FUNCTION
%   f: R^n+1 -> R; where the n+1 element is 't'. 
objfun = @(x) x(domdim + 1);

w1 = 0:0.01:1;
w2 = 1:0.01:2.5;
[W1, W2] = meshgrid(w1, w2);

x0 = zeros(1, domdim + 1);      % initial guess, the last element is 't'
ntheta = 2;                     % number of semi-infinite constraints
A = [];                         % equality constraints
b = [];                         % equality constraints 
Aeq = [];                       % inequality constraints
beq = [];                       % inequality constraints
lb = [0, 1];                    % 0 <= w1 and 0 <= w2
ub = [1, 2.5];                  % w1 <= 1 and w2 <= 2.5
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
    d = targetfun(W1, W2);

    % F(x, w)
    f = F(x(1:domdim), W1, W2);
    
    % t
    [cols, rows] = size(W1);
    t = zeros(cols, rows) + x(domdim + 1);

    % g_1(x, w)
    K1 = d - f - t;
    % g_2(x, w)
    K2 = f - d - t;
end

function apprxval = F(x, W1, W2)
    [c, r] = size(W1);
    apprxval = zeros(c, r);
    for i = 1:length(x)
        p1 = pows(1, i);
        p2 = pows(2, i);

        if p1 ~= 0 && p2 ~= 0
            apprxval = apprxval + x(i).*(W1.^pows(1, i)).*(W2.^pows(2, i));
        elseif p1 == 0
            apprxval = apprxval + x(i).*(W2.^pows(2, i));
        elseif p2 == 0
            apprxval = apprxval + x(i).*(W1.^pows(1, i));
        else 
            apprxval = apprxval + x(i);
        end
    end

end

function pows = getpows(d)
    if (d == 0)
        pows = zeros(2, 1);
        return;
    end
    pows = [zeros(2, d + 1), getpows(d - 1)];

    half = ceil(d/2);
    for i = 0:(half-1)
        % add (x^a)*(x^b)
        pows(1, i + 1) = d - i;
        pows(2, i + 1) = i;
        % add (x^b)*(x^a)
        pows(1, d - i + 1) = i;
        pows(2, d - i + 1) = d - i;
    end
    if rem(d, 2) == 0
        pows(1, half + 1) = half; 
        pows(2, half + 1) = half;
    end
end

end