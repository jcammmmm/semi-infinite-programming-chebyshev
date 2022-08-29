function capsiptest3d(targetfunlvl, degree, coeffs, domain)
    [W1, W2] = meshgrid(domain{1}, domain{2});

    pows = pnomialexps(3, degree);
    [c, r, d] = size(W1);
    approxlvl = @(W1, W2, lvl) pnomial3d(coeffs(1:length(pows)), {W1, W2, zeros(c, r, d) + lvl}, pows);

    % each time the plot is rotated 30 degrees around z axis
    for j = 0:60:360 
        for i = domain{3}
            cla
            target = targetfunlvl(W1, W2, i);
            approx = approxlvl(W1, W2, i);
            
            % ERROR PLOT
            % figure that contains the diference between target and the approx.
            tiledlayout(1, 2);
            graph1 = nexttile;
            imagesc(abs(target - approx));
            colormap(graph1, copper);
            colorbar;
            
            % APPROXIMATION PLOT
            graph2 = nexttile;
            % plot level sets target
            surf(W1, W2, target);
            colormap(graph2, "parula")
            hold on
            % plot level sets approximations
            mesh(W1, W2, approx);
            % view([45 + j, 45]);
            pause(0.2);
        end

    end

    % surf(W1, W2, W3); 
end