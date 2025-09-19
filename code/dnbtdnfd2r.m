%% 2 rotating vortex dynamics on non equispaced grid,  

clearvars;
close all;

%% space discretization

L = 20;

xc = 3;
yc = 0;
delta = 0.025;
h1 = 0.05;

[x, hx, mx2] = neqdeltagrid0(xc, L, delta, h1);
[y, hy, my2] = neqdeltagrid0(xc, L, delta, h1);

% x axis
x = unique([-flip(x), x]); x = x(:);
hx = [flip(hx), hx]; hx = hx(:);
mx = length(x);

y = unique([-flip(y), y]); y = y(:);
hy = [flip(hy), hy]; hy = hy(:);
my = length(y);

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

rho2vortex0 = rho(X, Y, xc, yc).*rho(X, Y, - xc, yc);
phase2vortex0 = atan2(Y - yc, X - xc) + atan2(Y - yc, X + xc);

U0 = sqrt(rho2vortex0).*exp(1i*phase2vortex0);

%% plot phase
figure(1);
lim = 15;
L = [-lim, lim];
surf(X,Y,angle(U0), 'EdgeColor','none');
colorbar; colormap('jet');
xlim(L); ylim(L);
title('Initial condition');
view(2);

%% Strang time splitting integration
tau = 0.1;
nt = 36/tau;

Ex = expm(1i*tau/4*Dxx);
Ey = expm(1i*tau/4*Dyy);

t = 0;
U = U0;
L = [-15,15];
k = 0;

tic;

while t < 36

    k = k + 1;

    U = Ex*U*Ey.';
    U = exp(1i*tau/2*(1 - abs(U).^2)).*U;
    U = Ex*U*Ey.';
    t = t + tau;

    if mod(k,round(nt/10))==0 % Plot ogni 1/10 dei passi totali
        fprintf('step %d / %d, max|U-U_prev| = %.3e\n', k, nt-1, max(abs(U(:) - U0(:))));
        
        figure(3); % Usa una nuova figura per l'animazione o sovrapposizioni
        subplot(1,2,1);
        surf(X,Y,abs(U), 'EdgeColor','none');
        colorbar; colormap('jet');
        xlim(L); ylim(L);
        title(sprintf('Ampiezza |U| a t = %.2f', k*tau));
        view(2);
        
        subplot(1,2,2);
        surf(X,Y,angle(U), 'EdgeColor','none');
        colorbar; colormap('jet');
        xlim(L); ylim(L);
        title(sprintf('Fase angle(U) a t = %.2f', k*tau));
        view(2);
        
        drawnow; % Aggiorna la figura per mostrare l'animazione
    end

end

toc;

%% plot U as a 3D surface
% Generate a 3D surface plot of U0
Ld = 20;
figure;
subplot(1,2,1);
surf(X, Y, abs(U0), 'EdgeColor', 'none');
colorbar; colormap('jet');
xlim([-Ld,Ld]); ylim([-Ld,Ld]);
title('Magnitude of \psi_0');
view(2);

% surface for U(tstar)
subplot(1,2,2);
surf(X, Y, abs(U), 'EdgeColor', 'none');
colorbar; colormap('jet');
xlim([-Ld,Ld]); ylim([-Ld,Ld]);
title('Magnitude of \psi_{22}');
view(2);

%% grid plot
% Definisci i nodi della griglia

figure;
hold on;

% Linee verticali (costanti in x)
plot(X', Y', 'k-', 'LineWidth', 0.1);

% Linee orizzontali (costanti in y)
plot(X, Y, 'k-', 'LineWidth', 0.1);

xlabel('x'); ylabel('y');
title('Rotating vortex grid');
xlim([X(1), X(end)]);
ylim([Y(1), Y(end)]);
axis equal;
hold off;


%% save
save('data/rotating_nfd.mat', 'U', 'X', 'Y');