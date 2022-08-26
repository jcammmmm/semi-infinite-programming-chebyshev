% SIP OBJECTIVE FUNCTION
objfun = @(x) x(7);
w1 = 0:0.1:1;
w2 = 1:0.1:2.5;
[W1, W2] = meshgrid(w1, w2);

% d(w): target function
d = log(W1 + W2).*sin(W1);
% surf(W1, W2, d); 

x0 = zeros(1, 11);     % initial guess
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

function apprxval = F1(x, W1, W2) 
    % F(x, w; 2): approximation function
    apprxval =          ...
        x(1) +          ...
        x(2).*W2 +      ...
        x(3).*W2.^2 +   ...
        x(4).*W1 +      ...
        x(5).*W1.*W2 +  ...
        x(6).*W1.^2;
end

function apprxval = F2(x, W1, W2) 
    % F(x, w; 2): approximation function
    apprxval =              ...
        x(1) +              ...
        x(2).*W2    +       ...
        x(3).*(W2.^2) +       ...
        x(4).*(W2.^3) +       ...
        x(5).*W1    +       ...
        x(6).*(W1.^2) +       ...
        x(7).*(W1.^3) +       ...
        x(8).*W1.*W2 +      ...
        x(9).*(W1.^2).*W2 +    ...
        x(10).*W1.*(W2.^2) ; 
end

% SEMI-INFINITE CONSTRAINTS DEFINITION
% Function signature interface is provided by fseminf authors.
% X := (x, t)
% S := an optional sampling interval
function [c, ceq, K1, K2, S] = seminfcon(x, S)
    % non-linear constraints
    c = [];              
    ceq = [];

    % semi-infinite constraint parameter
    w1 = 0:0.005:1;
    w2 = 1:0.005:2.5;
    [W1, W2] = meshgrid(w1, w2);

    % d(w): target function
    d = log(W1 + W2).*sin(W1);

    % F(x, w)
    f = F2(x(1:10), W1, W2);

    % t
    [cols, rows] = size(W1);
    t = zeros(cols, rows) + x(7);

    % g_1(x, w)
    K1 = d - f - t;
    % g_2(x, w)
    K2 = f - d - t;
end