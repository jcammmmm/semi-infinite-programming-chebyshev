clear all
f = @f2;
r = f(5);
r

[x1, x2] = ndgrid(0:0.1:1);
x = {x1, x2};
for i = x
    disp(i{1}(8, 9))
end

function y = f1(x) 
    y = x*x;
end

function y = f2(x)
    y = x*x*x;
end