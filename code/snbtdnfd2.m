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

%% initial condition
rhoVortex = rho(X, Y, 0,0);
phase = atan2(Y,X);

U0 = sqrt(rhoVortex).*exp(1i*phase);
U = U0;

%% time integration
tstar = 20;
nt = 450;
tau = tstar/(nt - 1);

Ex = expm(1i*tau/4*Dxx);
Ey = expm(1i*tau/4*Dyy);

for k = 1:nt - 1
    
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
%subplot(1,2,1);
surf(X, Y, abs(U0), 'EdgeColor', 'none');
colorbar; colormap(color);
xlabel('X-axis');
ylabel('Y-axis');
xlim([-L,L]); ylim([-L,L]);
zlabel('Magnitude of \psi_0');
axis equal tight;
view(2);

% surface for U(tstar)
surf(X, Y, abs(U), 'EdgeColor', 'none');
colorbar, colormap(color);
xlabel('X-axis');
ylabel('Y-axis');
xlim([-L,L]); ylim([-L,L]);
zlabel('Magnitude of \psi');
axis equal tight;
view(2);

vmin = min([abs(U(:)); abs(U0(:))]);
vmax = max([abs(U(:)); abs(U0(:))]);

subplot(1,2,1); clim([vmin vmax]);
subplot(1,2,2); clim([vmin vmax]);

%% save data
U_nfd = U;
X_nfd = X;
Y_nfd = Y;
save('nfd.mat', 'U_nfd', 'X_nfd', 'Y_nfd');

