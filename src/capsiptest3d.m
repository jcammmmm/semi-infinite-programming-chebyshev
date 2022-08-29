function capsiptest3d(targetfunlvl, degree, coeffs, domain)
    [W1, W2] = meshgrid(domain{1}, domain{2});

    pows = pnomialexps(3, degree);
    [c, r, d] = size(W1);
    approxlvl = @(W1, W2, lvl) pnomial3d(coeffs(1:length(pows)), {W1, W2, zeros(c, r, d) + lvl}, pows);

    for j = 0:30:180
        for i = domain{3}
            cla
            % plot level sets target
            surf(W1, W2, targetfunlvl(W1, W2, i));
            hold on
            % plot level sets approximations
            mesh(W1, W2, approxlvl(W1, W2, i));
            view([45 + j, 45]);
            pause(0.2);
    
        end

    end

    % surf(W1, W2, W3); 
end