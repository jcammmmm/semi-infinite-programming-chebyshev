example = 3;

% set figure's window position and size
set(gcf, 'Position',  [100, 100, 800, 350])

if example < 8
    switch example
        
        % :: 4.1 reemtsen
        case 1 % log(w1 + w2)*sin(w1)
            w1 = 0:0.05:1;
            w2 = 1:0.05:2.5;
            d = @(W) log(W{1} + W{2}).*sin(W{1});
        
        % :: 4.2 reemtsen
        case 2 % (1 + w1)^w2
            w1 = 0:0.05:1;
            w2 = 1:0.05:2.5;
            d = @(W) (1 + W{1}).^W{2};
        
        % :: 4.5 reemtsen
        case 3 % e^(w1*w1 - w1*w2)
            w1 = -1:0.05:1;
            w2 = -1:0.05:1;
            d = @(W) exp(W{1}.^2 - W{1}.*W{2});

        % :: 4.4 reemtsen
        case 4 % 1/(w1 + 2*w2 + 4)
            w1 = -1:0.05:1;
            w2 = -1:0.05:1;
            d = @(W) 1./(W{1} + 2.*W{2} + 4);
    
        % :: 4.6 reemtsen
        case 5 % sqrt(w1 + 2*w2 + 4)
            w1 = -1:0.05:1;
            w2 = -1:0.05:1;
            d = @(W) sqrt(W{1} + 2.*W{2} + 4);

        % :: matlab common example
        case 6 %
            w1 = -1:0.05:1;
            w2 = -1:0.05:1;
            d = @(W) sin(W{1}.^2).*cos(W{2}.^2);

        % :: matlab common example
        case 7 % back
            w1 = -1:0.05:1;
            w2 = -1:0.05:1;
            d = @(W) W{1}.*exp(-W{1}.^2 - W{2}.^2);

    end
    degree = 9;
    coeffs = capsip(d, degree, {w1, w2});
    capsiptest2d(d, degree, coeffs, {w1, w2});
else 

    switch example
        % :: 4.3 reemtsen
        case 8 % cos (w3*(1 + w1)^w2)
            w1 = 0:0.05:1;
            w2 = 1:0.05:2;
            w3 = 0:0.05:1;
            d    = @(W) cos(W{3}.*(1 + W{1}).^W{2});
            dlvl = @(W1, W2, lvl) cos(lvl.*(1 + W1).^W2);

        % :: 4.8 reemtsen
        case 9 % | log((w1w2 + 1)/(w1 + 0.5)) | x2^((x3 + 1)/2)
            w1 = 0:0.05:1;
            w2 = 0:0.05:1;
            w3 = 0:0.05:1;
            d    = @(W) abs(log((W{1}.*W{2} + 1)./(W{1} + 0.5))).*(W{2}.^((W{3} + 1)/2));
            dlvl = @(W1, W2, lvl) abs(log((W1.*W2 + 1)./(W1 + 0.5))).*(W2.^((lvl + 1)/2));
            
        % :: our proposal
        case 10 % w1 + w2 + w3
            w1 = 0:0.05:1;
            w2 = 0:0.05:1;
            w3 = 0:0.05:1;
            d    = @(W) W{1} + W{2} + W{3};
            dlvl = @(W1, W2, lvl) W1 + W2 + lvl;

    end
    degree = 2
    coeffs = capsip(d, degree, {w1, w2, w3});
    capsiptest3d(dlvl, degree, coeffs, {w1, w2, w3});
end


