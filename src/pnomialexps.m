function pows = pnomialexps(dim, deg)
    if dim == 2
        pows = [zeros(1, 2)];
        for j = 1:deg
            half = ceil(j/2);
            dpow = zeros(j + 1, 2);
            for i = 0:(half-1)
                % add (x^a)*(x^b)
                dpow(i + 1, 1) = j - i;
                dpow(i + 1, 2) = i;
                % add (x^b)*(x^a)
                dpow(j - i + 1, 1) = i;
                dpow(j - i + 1, 2) = j - i;
            end
            if rem(j, 2) == 0
                dpow(half + 1, 1) = half; 
                dpow(half + 1, 2) = half;
            end
            pows = [pows; dpow];
        end
    elseif dim == 3
        pows = [];
        for i = 0:deg
            for j = 0:deg
                for k = 0:deg
                    if i + j + k <= deg
                        pows = [pows; [i j k]];
                    end
                end
            end
        end
    else
        error('Multivariate-polynomials with dimension greater than 3 are not supported yet.');
end