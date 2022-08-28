clear all
%{
x1 = 0:0.1:1;
x2 = 1:0.1:2;
x3 = 2:0.1:3;
x4 = 4:0.1:4;

[s1, s2, s3, s4] = ngrid(x1, x2, x3, x4)
%}

[A, B, C, D, E] = ndgrid(0:0.1:1);
R = sin(A).*sin(B) + cos(C).*cos(D) + E;
R(1, 2, 3, 4, 5)
sin(0).*sin(0.1) + cos(0.2).*cos(0.3) + 0.4;

