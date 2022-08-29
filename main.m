d = @(W1, W2) log(W1 + W2).*sin(W1);
d = @(W1, W2) (1 + W1).^W2;
d = @(W1, W2) W1.*exp(-W1.^2 - W2.^2);
d = @(W1, W2) sin(W1.^2).*cos(W2.^2);

x = capsip(d, 5);
capsip_test(d, 5, x);

function capsip_test(targetfun, degree, coeffs)
    w1 = 0:0.1:1;
    w2 = 1:0.1:2.5;
    [W1, W2] = meshgrid(w1, w2);

    d = targetfun(W1, W2);
    surf(W1, W2, d); 
    hold on;
    
    pows = getpows(degree);
    f = F(coeffs(1:length(pows)), W1, W2, pows);
    mesh(W1, W2, f);
end