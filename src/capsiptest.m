function capsiptest(targetfun, degree, coeffs, domain)
    domdim = length(domain);
    if domdim == 2
        [W1, W2] = meshgrid(domain{1}, domain{2});
        W = {W1, W2};
        pnomial = @pnomial2d;
    elseif domdim == 3
        [W1, W2, W3] = meshgrid(domain{1}, domain{2}, domain{3});
        W = {W1, W2, W3};
        pnomial = @pnomial3d;
    end
    % The surface is the target function
    target = targetfun(W);
    
    surf(W1, W2, target); 
    hold on;
    
    pows = pnomialexps(2, degree);
    approx = pnomial(coeffs(1:length(pows)), W, pows);
    
    % The mesh is the Chevyshev polynomial approximation
    mesh(W1, W2, approx);
end