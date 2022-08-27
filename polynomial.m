

[W1, W2] = meshgrid(0:0.4:1, 0:0.4:1);
r1 = F(ones(1, 10), W1, W2);
r2 = F2(ones(1, 10), W1, W2);

abs(r1 - r2)
r1(3, 1) - r2(3, 1)

function apprxval = F2(x, W1, W2) 
    % F(x, w; 2): approximation function
    apprxval =              ...
        x(1) +              ...
        x(2).*W2    +       ...
        x(3).*(W2.^2) +       ...
        x(4).*(W2.^3) +       ...
        x(5).*W1    +       ...
        x(6).*(W1.^2) +       ...
        x(7).*(W1.^3) +       ...
        x(8).*W1.*W2 +      ...
        x(9).*(W1.^2).*W2 +    ...
        x(10).*W1.*(W2.^2) ; 
end

function pows = getpows(d)
    if (d == 0)
        pows = zeros(2, 1);
        return;
    end
    pows = [zeros(2, d + 1), getpows(d - 1)];

    half = ceil(d/2);
    for i = 0:(half-1)
        % add (x^a)*(x^b)
        pows(1, i + 1) = d - i;
        pows(2, i + 1) = i;
        % add (x^b)*(x^a)
        pows(1, d - i + 1) = i;
        pows(2, d - i + 1) = d - i;
    end
    if rem(d, 2) == 0
        pows(1, half + 1) = half; 
        pows(2, half + 1) = half;
    end
end
  
function apprxval = F(x, W1, W2)
    pows = getpows(3);
    [c, r] = size(W1);
    apprxval = zeros(c, r);
    for i = 1:length(x)
        p1 = pows(1, i);
        p2 = pows(2, i);

        if p1 ~= 0 && p2 ~= 0
            apprxval = apprxval + x(i).*(W1.^pows(1, i)).*(W2.^pows(2, i));
        elseif p1 == 0
            apprxval = apprxval + x(i).*(W2.^pows(2, i));
        elseif p2 == 0
            apprxval = apprxval + x(i).*(W1.^pows(1, i));
        else 
            apprxval = apprxval + x(i);
        end
    end

end