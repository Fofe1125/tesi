clearvars;
close all;

%% space rescalation
% xl -> 2/pi*atan(xl/alphal)
alpha = 10;

% define transformation and inverse
transform = @(xl) (2/pi)*atan(xl/alpha);
antitransform = @(yl) alpha*tan((pi/2)*yl);

%% space discretization of (-1,1)^3 (3D now)

epsbd = 1e-10;
my1 = 157; 
my2 = 157;
my3 = 11; 

h1 = 2/(my1 - 1); 
h2 = 2/(my2 - 1);
h3 = 2/(my3 - 1);

y1 = linspace(-1 + epsbd,1 - epsbd,my1)';
y2 = linspace(-1 + epsbd,1 - epsbd,my2)';
y3 = linspace(-1 + epsbd,1 - epsbd,my3)'; 

y1p = 2./(pi*alpha*(tan(pi/2*y1).^2 + 1));
y2p = 2./(pi*alpha*(tan(pi/2*y2).^2 + 1));
y3p = 2./(pi*alpha*(tan(pi/2*y3).^2 + 1));

y1pp = -(4*tan(pi/2*y1))./(pi*alpha^2*(tan(pi/2*y1).^2 + 1).^2);
y2pp = -(4*tan(pi/2*y2))./(pi*alpha^2*(tan(pi/2*y2).^2 + 1).^2);
y3pp = -(4*tan(pi/2*y3))./(pi*alpha^2*(tan(pi/2*y3).^2 + 1).^2); 

[Y1,Y2,Y3] = ndgrid(y1,y2,y3); 

% Laplacian 
[Dx1, Dxx1] = buildMatrix(my1, h1); 
[Dx2, Dxx2] = buildMatrix(my2, h2);
[Dx3, Dxx3] = buildMatrix(my3, h3);

A1 = 1i/2*(spdiags(y1p.^2, 0, my1, my1)*Dxx1) + ...
    1i/2*(spdiags(y1pp, 0, my1, my1)*Dx1);
A2 = 1i/2*(spdiags(y2p.^2, 0, my2, my2)*Dxx2) + ...
    1i/2*(spdiags(y2pp, 0, my2, my2)*Dx2);
A3 = 1i/2*(spdiags(y3p.^2, 0, my3, my3)*Dxx3) + ...
    1i/2*(spdiags(y3pp, 0, my3, my3)*Dx3);

X1 = antitransform(Y1); 
X2 = antitransform(Y2);
X3 = antitransform(Y3);

%% Initial condition
xc1 = 3; yc1 = 0;
xc2 = -3; yc2 = 0;

theta1 = atan2(X2 - yc1, X1 - xc1);
theta2 = atan2(X2 - yc2, X1 - xc2);

r1 = (X1 - xc1).^2 + (X2 - yc1).^2;
r2 = (X1 - xc2).^2 + (X2 - yc2).^2;

rho1 = rho3d(r1);
rho2 = rho3d(r2);

U0 = sqrt(rho1.*rho2).*exp(1i*(theta1 + theta2));
U = U0;

%% plot vortex
c = 0.2;

figure(1);
hold on;

p = patch(isosurface(X1,X2,X3,abs(U).^2,c));
% isonormals(X1,X2,X3,abs(U).^2,p);
set(p, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.6);

camlight('headlight');
axis equal
xlabel('x'); ylabel('y'); zlabel('z');

l = 30; L = [-l, l];
xlim(L); ylim(L); zlim(L)

view(3);
hold off;

%% Strang splitting con plot evolutivo
tau = 0.02; 
tstar = 80; 
numSteps = tstar/tau; % number of time steps

E1 = expm(tau*A1);
E2 = expm(tau*A2);
E3 = expm(tau*A3);

% figure;
% c = 0.2;
% p = patch(isosurface(X1,X2,X3,abs(U).^2,c));
% set(p, 'FaceColor', 'magenta', 'EdgeColor', 'none', 'FaceAlpha', 0.6);
% camlight('headlight'); lighting gouraud;
% l = 30; L = [-l, l];
% xlim(L); ylim(L); zlim(L);
% view(3);

figure;
hSurf = surf(X1(:,:,5),X2(:,:,5),abs(U(:,:,5)), 'EdgeColor','none');
colormap('jet'), colorbar;
l = 30; L = [-l, l];
xlim(L); ylim(L); zlim(L);
view(2);

tic;
for k = 1:numSteps
    % Strang splitting
    U = ndcovlt(U,{E1,E2,E3});
    U = exp(tau*(1i/2)*(1 - abs(U).^2)).*U;
    U = ndcovlt(U,{E1,E2,E3});

    % Aggiorna isosurface ogni N step (per velocità)
    if mod(k,10)==0

        % delete(p);
        % p = patch(isosurface(X1,X2,X3,abs(U).^2,c));
        % set(p, 'FaceColor', 'magenta', 'EdgeColor', 'none', 'FaceAlpha', 0.6);

        set(hSurf, 'ZData', abs(U(:,:,5)));

        title(['Time = ', num2str(k*tau, '%.2f')]);
        drawnow;
        
    end
end
toc;

% soluzione finale
U = U1 .* U2;

%% plot vortex
figure;
hold on;

p1 = patch(isosurface(X1,X2,X3,abs(U1).^2,c));
p2 = patch(isosurface(X1,X2,X3,abs(U2).^2,c));

isonormals(X1,X2,X3,abs(U2).^2,p1); isonormals(X1,X2,X3,abs(U2).^2,p2);

set(p1, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.6);
set(p2, 'Facecolor', 'Red', 'Edgecolor', 'none', 'FaceAlpha', 0.6);

camlight('headlight');
axis equal
xlabel('x'); ylabel('y'); zlabel('z');

l = 30; L = [-l, l];
xlim(L); ylim(L); zlim(L)

view(3);
hold off;

%% section