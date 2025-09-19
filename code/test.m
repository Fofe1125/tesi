% dynamic of two vortex 4th order finite difference

clearvars;
close all;

epsbd = 1e-8;
L = 20;
mx = 301; my = 301;
hx = 2*L/(mx - 1); hy = 2*L/(my - 1);
x = linspace(-L + epsbd,L - epsbd,mx);
y = linspace(-L + epsbd,L - epsbd,my);

[X,Y] = ndgrid(x,y);

% laplacian in mapped space
[~, Dxx] = buildMatrix(mx, hx);
[~, Dyy] = buildMatrix(my, hy);
Dxx(1,1:2) = [-2,2]/(hx^2);
Dxx(mx,mx-1:mx) = [2,-2]/(hx^2);
Dyy(1,1:2) = [-2,2]/(hy^2);
Dyy(my,my-1:my) = [2,-2]/(hy^2);

%% vortex

xc = 1;
yc = 0;

rhoVortex = rho(X, Y, xc,yc).*rho(X,Y,-xc,yc);
phaseVortex = atan2(Y - yc, X - xc) + atan2(Y-yc, X-xc);

U0 = sqrt(rhoVortex).*exp(1i*phaseVortex);

%% plot initial condition
figure(1);
lim = 15;
L = [-lim, lim];
surf(X,Y,abs(U0), 'EdgeColor','none');
colorbar; colormap('jet');
xlim(L); ylim(L);
title('Initial condition');
view(2);

%% Strang time integration
tstar = 20;
U = U0;
nt = 450;
tau = tstar/(nt - 1);

Ex = expm(1i*tau/4*Dxx);
Ey = expm(1i*tau/4*Dyy);

for k = 1:nt - 1
    
    U = Ex*U*Ey.';
    U = exp(tau*(1i/2)*(1 - abs(U).^2)).*U;
    U = Ex*U*Ey.';

end

%% plot final solution
figure(2);
lim = 15;
L = [-lim, lim];
surf(X,Y,abs(U), 'EdgeColor','none');
colorbar; colormap('jet');
xlim(L); ylim(L);
title('Initial condition');
view(2);