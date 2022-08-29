function apprxval = pnomial2d(x, W, pows)
% Evaluates a polynomial by doing a dot product.
%   x := array of coefficients.
%   W := cell array of 2d-dimensional meshes
%   pows := 2*n-1 matrix of 2d polynomial powers

    [c, r] = size(W{1});
    apprxval = zeros(c, r);
    for i = 1:length(x)
        p1 = pows(1, i);
        p2 = pows(2, i);

        if p1 ~= 0 && p2 ~= 0
            apprxval = apprxval + x(i).*(W{1}.^pows(1, i)).*(W{2}.^pows(2, i));
        elseif p1 == 0
            apprxval = apprxval + x(i).*(W{2}.^pows(2, i));
        elseif p2 == 0
            apprxval = apprxval + x(i).*(W{1}.^pows(1, i));
        else 
            apprxval = apprxval + x(i);
        end
    end

end
