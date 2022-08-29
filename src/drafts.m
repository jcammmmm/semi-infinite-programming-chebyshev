
d = @(W1, W2, level) cos(W1.*(1 + W2).^level);
for j = 0:30:180
    for i = 0:0.1:1
        [W1, W2] = ndgrid(0:0.01:1);
        mesh(W1, W2, d(W1, W2, i));
        view([45 + j, 45]);
        pause(0.1);
    end
end
