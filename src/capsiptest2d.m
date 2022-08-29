function capsiptest2d(targetfun, degree, coeffs, domain)
    [W1, W2] = meshgrid(domain{1}, domain{2});
    W = {W1, W2};
    pnomial = @pnomial2d;
    
    % The surface is the target function
    target = targetfun(W);
    surf(W1, W2, target); 
    hold on;
    pows = pnomialexps(2, degree);
    approx = pnomial(coeffs(1:length(pows)), W, pows);
    % The mesh is the Chevyshev polynomial approximation
    mesh(W1, W2, approx);
end