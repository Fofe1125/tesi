%% 2 marching vortex dynamics on non equispaced grid

clearvars;
close all;

%% space discretization

Lx = 10;
Ly = 20;

xc = 5;
yc = 10;
delta = 0.05;
h1 = 0.05;

x = neqdeltagrid0(xc,Lx, delta, h1);

x = unique([-flip(x), x]); x = x(:);
hx = diff(x);
mx = length(x);

my = 151;
y = linspace(-Ly,Ly,my)';
hy = diff(y);

[X, Y] = ndgrid(x, y);

%% matrix
dux = 2./(hx(1:mx-2) .* (hx(1:mx-2) + hx(2:mx-1))); % upper diag 
dmx = -2./(hx(1:mx-2).*hx(2:mx-1)); % main diag
dlx = 2./(hx(2:mx-1).*(hx(1:mx-2) + hx(2:mx-1))); % lower diag

duy = 2./(hy(1:my-2) .* (hy(1:my-2) + hy(2:my-1))); % upper diag 
dmy = -2./(hy(1:my-2).*hy(2:my-1)); % main diag
dly = 2./(hy(2:my-1).*(hy(1:my-2) + hy(2:my-1))); % lower diag

Dxx = spdiags([[dux;0;0], [0; dmx; 0], [0;0;dlx]], -1:1, mx,mx);
Dyy = spdiags([[duy;0;0], [0; dmy; 0], [0;0;dly]], -1:1, my,my);

% Homogeneus Neumann conditions 
Dxx(1,1:2) = [-2,2]/(hx(1)^2);
Dyy(1,1:2) = [-2,2]/(hy(1)^2);
Dxx(mx,mx-1:mx) = [2,-2]/(hx(mx-1)^2);
Dyy(my,my-1:my) = [2,-2]/(hy(my-1)^2);

%% initial condition
rho2vortex = rho(X, Y, xc, yc).*rho(X, Y, - xc, yc);
phase2vortex = atan2(Y - yc, X - xc) - atan2(Y - yc, X + xc);

U0 = sqrt(rho2vortex).*exp(1i*phase2vortex);

%% Strang time splitting integration
tstar = 30;
nt = 1850; 
tau = tstar/nt;

Ex = expm(1i*tau/4*Dxx);
Ey = expm(1i*tau/4*Dyy);

U = U0;

for i = 1:nt - 1

    U = Ex*U*Ey.';
    U = exp(1i*tau/2*(1 - abs(U).^2)).*U;
    U = Ex*U*Ey.';

    if i == 617
        U_10 = U;
    end
 
end
%% plot U as a 3D surface
% Generate a 3D surface plot of U0
figure;
subplot(1,3,1);
surf(X, Y, abs(U0), 'EdgeColor', 'none');
colorbar; colormap('jet');
title('Magnitude of \psi_0');
xlim([-Lx,Lx]); ylim([-Ly,Ly]);
view(2);

subplot(1,3,2);
surf(X, Y, abs(U_10), 'EdgeColor', 'none');
colorbar; colormap('jet');
title('Magnitude of \psi_{20}');
xlim([-Lx,Lx]); ylim([-Ly,Ly]);
view(2);

% surface for U(tstar)
subplot(1,3,3);
surf(X, Y, abs(U), 'EdgeColor', 'none');
colorbar; colormap('jet');
title('Magnitude of \psi_{30}');
xlim([-Lx,Lx]); ylim([-Ly,Ly]);
view(2);

%% save matrix
save('data/straight_30t5.mat', 'U', 'X', 'Y');

%% grid plot
% Definisci i nodi della griglia

figure;
hold on;

% Linee verticali (costanti in x)
plot(X', Y', 'k-', 'LineWidth', 0.1);

% Linee orizzontali (costanti in y)
plot(X, Y, 'k-', 'LineWidth', 0.1);

xlabel('x'); ylabel('y');
title('Marching vortex grid');
xlim([X(1), X(end)]);
ylim([Y(1), Y(end)]);
axis equal;
hold off;