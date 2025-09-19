clearvars;
close all;

L = 20;
mx = 301;
my = 301;

x = linspace(-L,L,mx)'; hx = x(2) - x(1);
y = linspace(-L,L,my)'; hy = y(2) - y(1);

[~, Dxx] = buildMatrix(mx,hx);
[~, Dyy] = buildMatrix(my,hy);

