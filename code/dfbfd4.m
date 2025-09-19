% dynamic of two vortex 4th order finite difference

clearvars;
close all;

%% space rescalation
% xl -> 2/pi*atan(xl/alphal)
alpha = 10;

% define transformation and inverse
transform = @(xl) (2/pi)*atan(xl/alpha);
antitransform = @(yl) alpha*tan((pi/2)*yl);

%% space discretization of (-1,1)^2

epsbd = 1e-10;
my1 = 101; my2 = 101;
h = 2/(my1 - 1);
y1 = linspace(-1 + epsbd,1 - epsbd,my1)';
y2 = linspace(-1 + epsbd,1 - epsbd,my2)';

x1 = antitransform(y1);
x2 = antitransform(y2);

% y1p = 2*alpha./(pi*(x1).^2 + alpha^2);
% y2p = 2*alpha./(pi*(x2).^2 + alpha^2);
% y1pp = -(4*alpha*y1)./(pi*(x1.^2+alpha.^2).^2);
% y2pp = -(4*alpha*y2)./(pi*(x2.^2+alpha.^2).^2);
y1p = 2./(pi*alpha*(tan(pi/2*y1).^2 + 1));
y2p = 2./(pi*alpha*(tan(pi/2*y2).^2 + 1));
y1pp = -(4*tan(pi/2*y1))./(pi*alpha^2*(tan(pi/2*y1).^2 + 1).^2);
y2pp = -(4*tan(pi/2*y2))./(pi*alpha^2*(tan(pi/2*y2).^2 + 1).^2);

[Y1,Y2] = ndgrid(y1,y2);

% laplacian in mapped space
[Dx, Dxx] = buildMatrix(my1, h);
A1 = 1i/2*(spdiags(y1p.^2, 0, my1, my1)*Dxx) + ...
    1i/2*(spdiags(y1pp, 0, my1, my1)*Dx);
A2 = 1i/2*(spdiags(y2p.^2, 0, my2, my2)*Dxx) + ...
    1i/2*(spdiags(y2pp, 0, my2, my2)*Dx);

%% vortex 
X1 = antitransform(Y1); X2 = antitransform(Y2);

xc = 3;
yc = 0;

rho2vortex0 = rho(X1, X2, xc, yc).*rho(X1,X2,-xc,-yc);
phase2vortex0 = atan2(X2 - yc, X1 - xc) + atan2(X2 + yc, X1 + xc);

U0 = sqrt(rho2vortex0).*exp(1i*phase2vortex0);
%% Strang time integration

U = U0;

% tstar = 2;
% nt = 450;
% tau = tstar/(nt - 1);
tau = 0.1;
tstar = 80;
t = 0;

E1 = expm(tau*A1);
E2 = expm(tau*A2);

% for k = 1:nt - 1
while t < tstar
    
    U = E1*U*E2.';
    U = exp(tau*(1i/2)*(1 - abs(U).^2)).*U;
    U = E1*U*E2.';
    t = t + tau;
 
end

%% plot final solution
figure(2);

load('stepper_cmap.mat', 'CustomColormap');

lim = 10;
L = [-lim, lim];

subplot(1,2,1);
surf(X1,X2,abs(U), 'EdgeColor','none');
colorbar; colormap(CustomColormap);
xlim(L); ylim(L);
title('sol');
view(2);

subplot(1,2,2);
surf(X1,X2,abs(U0), 'EdgeColor','none');
colorbar; colormap('jet');
xlim(L); ylim(L);
title('U0');
view(2);

vmin = min([abs(U(:)); abs(U0(:))]);
vmax = max([abs(U(:)); abs(U0(:))]);

subplot(1,2,1); clim([vmin vmax]);
subplot(1,2,2); clim([vmin vmax]);
%% grid plot
% Definisci i nodi della griglia

% figure;
% hold on;
% 
% lim = 50;
% L = [-lim, lim];
% 
% % Linee verticali (costanti in x)
% plot(X1', X2', 'k-', 'LineWidth', 0.1);
% 
% % Linee orizzontali (costanti in y)
% plot(X1, X2, 'k-', 'LineWidth', 0.1);
% 
% xlabel('x'); ylabel('y');
% title('Rotating vortex grid');
% xlim(L);
% ylim(L);
% axis equal;
% hold off;