function apprxval = pnomial3d(x, W, pows)
% Evaluates a polynomial by doing a dot product.
%   x := array of coefficients.
%   W := cell array of 2d-dimensional meshes
%   pows := 2*n-1 matrix of 2d polynomial powers

    [c, r] = size(W{1});
    apprxval = zeros(c, r);
    for i = 1:length(x)
        p1 = pows(1, i);
        p2 = pows(2, i);
        p3 = pows(3, i);
        apprxval = apprxval + x(i).*(W{1}.^pows(1, i)).*(W{2}.^pows(2, i));
    end

end
