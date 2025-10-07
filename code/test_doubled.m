% Stationary Neumann boundary truncated domain finite difference 2nd order
% non uniform grid

%% non uniform grid
clearvars;
close all;

tic;
L = 30;
delta = 0.02;
h1 = 0.02;

[xhalf, ~] = neqdeltagrid0(0,L,delta,h1);
[yhalf, ~] = neqdeltagrid0(0,L,delta,h1);

x = unique([-flip(xhalf), xhalf]);
y = unique([-flip(yhalf), yhalf]);

hx = diff(x);
hy = diff(y);

x = x(:); y = y(:);
hx = hx(:); hy = hy(:);

mx = length(x);
my = length(y);

[X, Y] = ndgrid(x,y);

% matrixes
duxx = 2./(hx(1:mx-2) .* (hx(1:mx-2) + hx(2:mx-1))); % upper diag 
dmxx = -2./(hx(1:mx-2).*hx(2:mx-1)); % main diag
dlxx = 2./(hx(2:mx-1).*(hx(1:mx-2) + hx(2:mx-1))); % lower diag

duyy = 2./(hy(1:my-2) .* (hy(1:my-2) + hy(2:my-1))); 
dmyy = -2./(hy(1:my-2).*hy(2:my-1)); 
dlyy = 2./(hy(2:my-1).*(hy(1:my-2) + hy(2:my-1)));

Dxx = spdiags([[duxx;0;0], [0;dmxx;0], [0;0;dlxx]], -1:1, mx,mx);
Dyy = spdiags([[duyy;0;0], [0;dmyy;0], [0;0;dlyy]], -1:1, my, my);

% Homogeneus Nuemann BC
Dxx(1,1:2) = [-2 2]/(hx(1)^2); Dxx(mx,mx-1:mx) = [2 -2]/(hx(mx-1)^2);
Dyy(1,1:2) = [-2 2]/(hy(1)^2); Dyy(my,my-1:my) = [2 -2]/(hy(my-1)^2);

%% save grid
save('data/grid_nfd.mat', 'X', 'Y');

%% initial condition
rhoVortex = rho(X, Y, 0,0);
phase = atan2(Y,X);

U0 = sqrt(rhoVortex).*exp(1i*phase);
U = U0;

%% time integration
tstar = 20;
nt = 450;
tau = tstar/(nt - 1);
disp(['tau = ', num2str(tau)]);

Ex = expm(1i*tau/4*Dxx);
Ey = expm(1i*tau/4*Dyy);

tol = 1e-8;
mask = abs(U0) > tol;
counter = 1;

E = NaN(floor(nt/10),1);
M = NaN(floor(nt/10),1);

tic;
for k = 1:nt
    
    U = Ex*U*Ey.';
    U = exp(tau*(1i/2)*(1 - abs(U).^2)).*U;
    U = Ex*U*Ey.';
 
end

toc;

%% plot 
load('stepper_cmap.mat', 'CustomColormap');
color = CustomColormap;

L = 30;
figure;
subplot(1,2,1);
surf(X, Y, abs(U0), 'EdgeColor', 'none');
colorbar; colormap(color);
xlabel('X-axis');
ylabel('Y-axis');
xlim([-L,L]); ylim([-L,L]);
zlabel('Magnitude of \psi_0');
view(2);

subplot(1,2,2);
surf(X, Y, abs(U), 'EdgeColor', 'none');
colorbar, colormap(color);
xlabel('X-axis');
ylabel('Y-axis');
xlim([-L,L]); ylim([-L,L]);
zlabel('Magnitude of \psi');
view(2);

vmin = min([abs(U(:)); abs(U0(:))]);
vmax = max([abs(U(:)); abs(U0(:))]);

subplot(1,2,1); clim([vmin vmax]);
subplot(1,2,2); clim([vmin vmax]);

