clear all;

a = -1;
b = 1;
tol= 1e-7;
nrif=200;%100;
n= nrif;
sx = a;
dx = b;
xrif = linspace(sx,dx,n+2)';
x0 = xrif(2:n+1);
h = (dx-sx)/(n+1);
A = -spdiags(ones(n,1)*[1,-2,1]/h^2,-1:1,n,n);
b = ones(n,1);
F = @(u) 1/2*u'*A*u-b'*u; 
JF = @(u)A*u-b;
[x, its] = conjugate_gradient(x0,F,JF,sx,dx,tol);
err_rel = norm(A\b-x)/norm(A\b);