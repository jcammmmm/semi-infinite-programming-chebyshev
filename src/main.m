example = 5;

if example < 5
    w1 = 0:0.1:1;
    w2 = 1:0.1:2.5;
    switch example

        case 1
            d = @(W) log(W{1} + W{2}).*sin(W{1});
        case 2
            d = @(W) (1 + W{1}).^W{2};
        case 3
            d = @(W) W{1}.*exp(-W{1}.^2 - W{2}.^2);
        case 4
            d = @(W) sin(W{1}.^2).*cos(W{2}.^2);

    end
    coeffs = capsip(d, 3, {w1, w2});
    capsiptest2d(d, 3, coeffs, {w1, w2});
else 
    w1 = 0:0.01:1;
    w2 = 0:0.01:1;
    w3 = 0:0.01:1;
    switch example
        case 5
            d = @(W) cos(W{3}.*(1 + W{1}).^W{2});
            dlvl = @(W1, W2, lvl) cos(lvl.*(1 + W1).^W2);
    end
    coeffs = capsip(d, 2, {w1, w2, w3});
    capsiptest3d(dlvl, 2, coeffs, {w1, w2, w3});
end


