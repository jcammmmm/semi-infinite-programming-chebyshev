example = 7;

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
    w1 = 0:0.05:1;
    w2 = 0:0.05:1;
    w3 = 0:0.05:1;
    switch example
        case 5 
            % cos (w3*(1 + w1)^w2)
            d    = @(W) cos(W{3}.*(1 + W{1}).^W{2});
            dlvl = @(W1, W2, lvl) cos(lvl.*(1 + W1).^W2);
        case 6
            % w1 + w2 + w3
            d    = @(W) W{1} + W{2} + W{3};
            dlvl = @(W1, W2, lvl) W1 + W2 + lvl;
        case 7
            % | log((w1w2 + 1)/(w1 + 0.5)) | x2^((x3 + 1)/2)
            d    = @(W) abs(log((W{1}.*W{2} + 1)./(W{1} + 0.5))).*(W{2}.^((W{3} + 1)/2));
            dlvl = @(W1, W2, lvl) abs(log((W1.*W2 + 1)./(W1 + 0.5))).*(W2.^((lvl + 1)/2));
    end
    coeffs = capsip(d, 2, {w1, w2, w3});
    capsiptest3d(dlvl, 2, coeffs, {w1, w2, w3});
end


