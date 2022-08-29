function apprxval = pnomial3d(x, W, pows)
% Evaluates a polynomial by doing a dot product.
%   x := array of coefficients.
%   W := cell array of 2d-dimensional meshes
%   pows := n-1 x 2 matrix of 3-dimensional polynomial powers

    [c, r] = size(W{1});
    apprxval = zeros(c, r);
    for i = 1:length(x)
        p1 = pows(i, 1);
        p2 = pows(i, 2);
        p3 = pows(i, 3);
        apprxval = apprxval + x(i).*(W{1}.^p1).*(W{2}.^p2).*(W{3}.^p3);
    end

end
