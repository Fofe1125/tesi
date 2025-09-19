% Neumann Bc truncated domain finite diff 2-nd order equispaced
% 2 vortex

clearvars;
close all;

tic;

% space discretization
Lx = -20; Rx = 20;
Ly = -20; Ry = 20;
mx = 771;
my = 771;
hx = (Rx - Lx)/(mx - 1);
hy = (Ry - Ly)/(my - 1);

Dxx = spdiags(ones(mx,1)*[1, -2, 1]/(hx^2), -1:1, mx,mx);
Dyy = spdiags(ones(my,1)*[1, -2, 1]/(hy^2), -1:1, my,my);

% Homogeneus Neumann conditions 
Dxx(1,1:2) = [-2,2]/(hx^2);
Dyy(1,1:2) = [-2,2]/(hy^2);
Dxx(mx,mx-1:mx) = [2,-2]/(hx^2);
Dyy(my,my-1:my) = [2,-2]/(hy^2);

epsbd = 1e-8;
x = linspace(Lx+epsbd, Rx-epsbd, mx)';
y = linspace(Ly+epsbd, Ry-epsbd, my)';

xc1 = -3; yc1 = 0;
xc2 = 3; yc2 = 0;

[X, Y] = ndgrid(x, y);

rho2vortex = rho(X, Y, xc1, yc1).*rho(X, Y, xc2, yc2);
phase2vortex = atan2(Y - yc1, X - xc1) + atan2(Y - yc2, X - xc2);

U0 = sqrt(rho2vortex).*exp(1i*phase2vortex);

%% Strang time splitting integration
tau = 0.05;
t = 0;

Ex = expm(1i*tau/4*Dxx);
Ey = expm(1i*tau/4*Dyy);

U = U0;

while t < 36

    U = Ex*U*Ey.';
    U = exp(1i*tau/2*(1 - abs(U).^2)).*U;
    U = Ex*U*Ey.';
    t = t + tau;
    
end
toc;

%% plot 
figure;
subplot(1,2,1);
surf(X, Y, abs(U0), 'EdgeColor', 'none');
colorbar; colormap('jet');
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Magnitude of \psi_0');
title('3D Surface Plot of Magnitude of \psi_0');
view(2);

% surface for U(tstar)
subplot(1,2,2);
surf(X, Y, abs(U), 'EdgeColor', 'none');
colorbar; colormap('jet');
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Magnitude of \psi');
title('3D Surface Plot of Magnitude of \psi at t = ', num2str(t));
view(2);

%% save
save('data/rotating_ufd.mat', 'U', 'X', 'Y');