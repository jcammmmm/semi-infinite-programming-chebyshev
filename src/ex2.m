% SIP OBJECTIVE FUNCTION
objfun = @(x) x(7);
w1 = 0:0.1:1;
w2 = 1:0.1:2.5;
[W1, W2] = meshgrid(w1, w2);

% d(w): target function
d = (1 + W1).^W2;
surf(W1, W2, d); 
hold on;

x = [0.0000    1.0000   -1.1706    0.3176   -1.8111    0.9539    0.1471    2.1799   -0.4525   -0.4352         0]
x = [0.3309    1.0689   -0.5097    0.0733    2.4961   -1.8813    0.0366   -2.0208    1.5337    0.8357         0]
x = [0.3651    1.0000   -0.4667    0.0650    2.5107   -1.9049    0.0536   -2.0308    1.5325    0.8389    0.0366]
x = [0.0506    1.5327    0.8398    0.0626   -1.9007   -2.0342   -0.4544    2.5117    0.9804    0.3749    0.0366]
x = [0.3651    1.0000   -0.4667    0.0650    2.5107   -1.9049    0.0536   -2.0308    1.5325    0.8389    0.0366] % F2
x = [0.3651    1.0000   -0.4667    0.0650    2.5107   -1.9049    0.0536   -2.0308    1.5325    0.8389    0.0366] % F2
x = [0.0506    1.5327    0.8398    0.0626   -1.9007   -2.0342   -0.4544    2.5117    0.9804    0.3749    0.0366] % F
x = [0.0034    1.5359    0.8516    0.0295   -1.8358   -2.0786   -0.2850    2.5236    0.7099    0.5087    0.0369]
x = [-0.3138   -0.2614   -0.0590    0.0111    1.1002    0.7818   -0.0691   -0.6646    0.1342   -0.0796    0.0035]
x = [   -0.0000    1.0000   -0.3201    0.0418   -1.5528    0.4341   -0.0850    1.0220    0.0484   -0.1240    0.1188]
F = Fn(x(1:10), W1, W2);
mesh(W1, W2, F);

function apprxval = Fn(x, W1, W2)
    disp('ok')
    pows = getpows(3);
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

