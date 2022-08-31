exps = [];
d = 2; 
for i = 0:d
    for j = 0:d
        for k = 0:d
            if i + j + k <= d
                exps = [exps; [i j k]];
            end
        end
    end
end
exps