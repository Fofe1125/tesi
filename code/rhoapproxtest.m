clear all;
close all;

%% finite differences
m = 5001;
s = linspace(eps,1-eps,m)';
h = s(2) - s(1);

D1 = spdiags(ones(m,1)*[-1 0 1]/(2*h), -1:1, m,m);
D2 = spdiags(ones(m,1)*[1 -2 1]/(h^2), -1:1, m,m);

D1(1,:) = 0; D1(m,:) = 0;
D2(1,1:2) = [1 0]/(h^2); D2(m,m-1:m) = [0 1]/(h^2);

fun = @(g) ((s-1).^4).*(D2*g) + (2*(s-1).^3).*(D1*g) - ...
    (((s-1).^3)./s).*(D1*g) - (((s-1).^2)./(s.^2)).*g + ...
    (1-g.^2).*g;

Jfun = @(g) spdiags((s-1).^4, 0, m,m)*D2 + spdiags(2*(s-1).^3, 0, m,m)*D1-...
    spdiags(((s-1).^3)./s, 0, m,m)*D1 - spdiags(((s-1).^2)./(s.^2), 0, m,m) + ...
    speye(m) - spdiags(3*g.^2, 0, m,m);

g0 = linspace(0,1,m)';
g = g0;

tol = 1e-8;
delta = -Jfun(g) \ fun(g);
its = 1;

while norm(delta, inf) > tol

    g = g + delta;
    delta = -Jfun(g) \ fun(g);

    its = its + 1; % Increment iteration count

end

disp(['Newton converged in ', num2str(its), ' iterations'])

g = g + delta;

% retrive rho
r = s./(1 - s);
ffd = g.^2;

%% pade 2,3 and 4
r1 = linspace(0,100,300);
f2 = rho2(r1);
f3 = rho3(r1);
f4 = rho(r1,0,0,0);

%% plot results
figure;

subplot(1,2,1);
hold on;
plot(r,ffd,'.');
xlim([0 10]);
legend('FD, m = 1001');
grid on;
hold off;

subplot(1,2,2);
hold on;
plot(r,ffd);
plot(r,ones(length(r),1),'--','Color',[0.62, 0.23, 0.96]);
plot(r1,f2,'gs');
xlim([8 25]);
grid on;



