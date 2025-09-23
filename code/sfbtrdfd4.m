% stationary vortex simulation with
% trasformated domain free-boundary finite difference 4th order
% time splitting 

clear all;
close all;

%% space rescalation
% xl -> 2/pi*atan(xl/alphal)
alpha = 5;

% define transformation and inverse
transform = @(xl) (2/pi)*atan(xl/alpha);
antitransform = @(yl) alpha*tan((pi/2)*yl);

%% space discretization of (-1,1)^2

epsbd = 1e-8;
my1 = 501; my2 = 501;
h = 2/(my1 - 1);
y1 = linspace(-1 + epsbd,1 - epsbd,my1)';
y2 = linspace(-1 + epsbd,1 - epsbd,my2)';

y1p = 2./(pi*alpha*(tan(pi/2*y1).^2 + 1));
y2p = 2./(pi*alpha*(tan(pi/2*y2).^2 + 1));
y1pp = -(4*tan(pi/2*y1))./(pi*alpha^2*(tan(pi/2*y1).^2 + 1).^2);
y2pp = -(4*tan(pi/2*y2))./(pi*alpha^2*(tan(pi/2*y2).^2 + 1).^2);

[Y1,Y2] = ndgrid(y1,y2);

% laplacian in mapped space
[Dx, Dxx] = buildMatrix(my1, h);
A1 = 1i/2*(spdiags(y1p.^2, 0, my1, my1))*Dxx + ...
    1i/2*(spdiags(y1pp, 0, my1, my1))*Dx;
A2 = 1i/2*(spdiags(y2p.^2, 0, my2, my2))*Dxx + ...
    1i/2*(spdiags(y2pp, 0, my2, my2))*Dx;

%% vortex 
X1 = antitransform(Y1); X2 = antitransform(Y2);
% rhoVortex = rho(Y1, Y2, 0, 0);
% phaseVortex = atan2(Y2, Y1);
rhoVortex = rho(X1, X2, 0,0);
phaseVortex = atan2(X2,X1);

U0 = sqrt(rhoVortex).*exp(1i*phaseVortex);

%% time splitting integration
tstar = 20;
U = U0;
nt = 100;
tau = tstar/(nt - 1);

E1 = expm(tau/2*A1);
E2 = expm(tau/2*A2);

for k = 1:nt - 1
    
    U = E1*U*E2.';
    U = exp(tau*(1i/2)*(1 - abs(U).^2)).*U;
    U = E2*U*E1.';
 
end

Uref = U;

%% plot 
L = 10;
figure;
subplot(1,2,1);
surf(X1, X2, abs(U0), 'EdgeColor', 'none');
colorbar; colormap('jet');
xlabel('X-axis');
ylabel('Y-axis');
xlim([-L,L]); ylim([-L,L]);
zlabel('Magnitude of \psi_0');
view(2);

% surface for U(tstar)
subplot(1,2,2);
surf(X1, X2, abs(U), 'EdgeColor', 'none');
colorbar; colormap('jet');
xlabel('X-axis');
ylabel('Y-axis');
xlim([-L,L]); ylim([-L,L]);
zlabel('Magnitude of \psi');
view(2);

