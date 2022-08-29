d = @(W) log(W{1} + W{2}).*sin(W{1});
d = @(W) (1 + W{1}).^W{2};
d = @(W) W{1}.*exp(-W{1}.^2 - W{2}.^2);
d = @(W) sin(W{1}.^2).*cos(W{2}.^2);
d = @(W) cos(W{3}.*(1 + W{1}).^W{2});

w1 = 0:0.01:1;
w2 = 1:0.01:2.5;
w3 = 0:0.01:1;

x = capsip(d, 3, {w1, w2});
capsip_test(d, 3, x);

function capsip_test(targetfun, degree, coeffs)
    w1 = 0:0.1:1;
    w2 = 1:0.1:2.5;
    [W1, W2] = meshgrid(w1, w2);
    W = {W1, W2};

    d = targetfun(W);
    surf(W1, W2, d); 
    hold on;
    
    pows = pnomialexps(degree);
    f = pnomial2d(coeffs(1:length(pows)), W, pows);
    mesh(W1, W2, f);
end