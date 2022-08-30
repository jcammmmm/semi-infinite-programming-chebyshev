function capsiptest2d(targetfun, degree, coeffs, domain)
    [W1, W2] = meshgrid(domain{1}, domain{2});
    W = {W1, W2};
    
    target = targetfun(W);
    pows = pnomialexps(2, degree);
    approx = pnomial2d(coeffs(1:length(pows)), W, pows);
    
    tiledlayout(1, 2);
    graph1 = nexttile;
    imagesc(abs(target - approx));
    colormap(graph1, copper);
    colorbar;

    graph2 = nexttile;
    for j = 0:1:105
        cla
        view (j, 45)
        % The surface is the target function
        surf(W1, W2, target); 
        colormap(graph2, "parula")
        hold on
        % The mesh is the Chevyshev polynomial approximation
        mesh(W1, W2, approx);
        
        pause(0.01);
    end


end