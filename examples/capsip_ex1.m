% SIP OBJECTIVE FUNCTION
objfun = @(x) x(7);
w1 = 0:0.1:1;
w2 = 1:0.1:2.5;
[W1, W2] = meshgrid(w1, w2);

% d(w): target function
d = log(W1 + W2).*sin(W1);
surf(W1, W2, d); 
hold on;

x = [ 0.0000    1.0000   -1.1706    0.3176   -1.8111    0.9539    0.1471    2.1799   -0.4525   -0.4352         0]
x = [ 0.0000    1.0000   -0.2414   -0.0052   -1.5528    0.1589    0.1188    1.2186   -0.2071   -0.0254         0]
x = [ 0.0000    1.0000   -0.3160    0.0428   -1.5528    0.4198   -0.0894    1.0322    0.0541   -0.1263    0.1188]
x = [-0.0000    1.0000   -0.3206    0.0417   -1.5528    0.4361   -0.0845    1.0206    0.0473   -0.1233    0.1188]
x = [-0.0000    1.0000   -0.3201    0.0418   -1.5528    0.4341   -0.0850    1.0220    0.0484   -0.1240    0.1188]
x = [-0.0000    1.0000   -0.3201    0.0418   -1.5528    0.4341   -0.0850    1.0220    0.0484   -0.1240    0.1188] % 
x = [-0.3141   -0.2613   -0.0590    0.0109    1.1005    0.7816   -0.0686   -0.6645    0.1333   -0.0791    0.0035] % modify lower bounds
x = [-0.0000    1.0000   -0.3201    0.0418   -1.5528    0.4341   -0.0850    1.0220    0.0484   -0.1240    0.1188] % with bad lower bounds
x = [0.0333    0.0000    0.0012    0.0000    0.1831    0.3500         0    0.0007    0.0000    0.0000    0.0454]
x = [-0.3141   -0.2613   -0.0590    0.0109    1.1005    0.7816   -0.0686   -0.6645    0.1333   -0.0791    0.0035] % fix again lower bounds
x = [-0.3140   -0.2613   -0.0590    0.0110    1.1005    0.7817   -0.0690   -0.6646    0.1339   -0.0794    0.0035] % added upper bounds
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