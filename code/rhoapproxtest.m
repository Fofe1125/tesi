% solve_compactified.m
clearvars; close all; clc

% discretizzazione su s in [0,1]
m = 5001;
s = linspace(0,1,m)';      % s(1)=0, s(end)=1
h = s(2)-s(1);

% derivata prima e seconda (matrici centrali)
e = ones(m,1);
D1 = spdiags(e*[-1,0,1]/(2*h), -1:1, m,m);
D2 = spdiags(e*[1,-2,1]/(h^2), -1:1, m,m);

D1(1,:) = 0; D1(m,:) = 0;
D2(1,1:2) = 0; D2(m,m-1:m) = 1;

nonlinear_term = @(g) (1 - g.^2).*g;
nonlinear_Jdiag = @(g) (1 - 3*g.^2);   
Ffun = @(g) ( (s-1).^4 ).* (D2*g) + 2*(s-1).^3 .* (D1*g) - (((s-1).^3)./s).* (D1*g) ...
           - (((s-1).^2)./(s.^2)).* g + nonlinear_term(g);

Jfun = @(g) spdiags((s-1).^4,0,m,m)*D2 + spdiags(2*(s-1).^3,0,m,m)*D1 ...
           - spdiags(((s-1).^3)./s,0,m,m)*D1 - spdiags(((s-1).^2)./(s.^2),0,m,m) ...
           + spdiags(nonlinear_Jdiag(g),0,m,m);

g = (linspace(0,1,m)').^2; 
g(1) = 0;
g(end) = 1;

tol = 1e-8;
maxit = 100;
for k = 1:maxit
    F = Ffun(g);
    J = Jfun(g);

    F(1) = g(1) - 0;
    F(end) = g(end) - 1;
    
    delta = - (J \ F);
    g = g + delta;

    if norm(delta, Inf) < tol
        fprintf('Newton converged in %d iters, ||delta||_inf = %.3e\n', k, norm(delta,Inf));
        break;
    end

end

%% interpolation
pp = spline(s,g);
x = linspace(9,25,150);

r = linspace(0,20,150);
rhofd = ppval(pp,r./(1 + r)).^2;

r2 = linspace(0,20,20);
r3 = linspace(1,20,20);
r4 = linspace(2,20,20);

figure(1);
hold on;
plot(r,rhofd,'k');
plot(r2,rho2(r2), 'bo');
plot(r4,rho(r4,0,0,0),'r^');

legend('FD, N = 5001', 'q = 2', 'q = 4', 'Location','southeast');
grid on;

hold off;

figure(2);
hold on;
plot(x,ones(length(x)), 'm--');
plot(x,rho2(x), 'bo');
plot(x,rho(x,0,0,0), 'r^');
hold off;
grid on;
legend('1','q = 2', 'q = 4', 'Location','southeast');


%% pade approximations

