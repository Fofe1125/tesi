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
my1 = 101; 
my2 = 101;
my3 = 101; 
h = 2/(my1 - 1); 

y1 = linspace(-1 + epsbd,1 - epsbd,my1)';
y2 = linspace(-1 + epsbd,1 - epsbd,my2)';
y3 = linspace(-1 + epsbd,1 - epsbd,my3)'; % New dimension

x1 = antitransform(y1);
x2 = antitransform(y2);
x3 = antitransform(y3); % New dimension

% Derivatives for mapped space (assuming they are independent per dimension)
y1p = 2./(pi*alpha*(tan(pi/2*y1).^2 + 1));
y2p = 2./(pi*alpha*(tan(pi/2*y2).^2 + 1));
y3p = 2./(pi*alpha*(tan(pi/2*y3).^2 + 1));
y1pp = -(4*tan(pi/2*y1))./(pi*alpha^2*(tan(pi/2*y1).^2 + 1).^2);
y2pp = -(4*tan(pi/2*y2))./(pi*alpha^2*(tan(pi/2*y2).^2 + 1).^2);
y3pp = -(4*tan(pi/2*y3))./(pi*alpha^2*(tan(pi/2*y3).^2 + 1).^2); 

[Y1,Y2,Y3] = ndgrid(y1,y2,y3); 

% Laplacian in mapped space 
[Dx, Dxx] = buildMatrix(my1, h); 

A1 = 1i/2*(spdiags(y1p.^2, 0, my1, my1)*Dxx) + ...
    1i/2*(spdiags(y1pp, 0, my1, my1)*Dx);
A2 = 1i/2*(spdiags(y2p.^2, 0, my2, my2)*Dxx) + ...
    1i/2*(spdiags(y2pp, 0, my2, my2)*Dx);
A3 = 1i/2*(spdiags(y3p.^2, 0, my3, my3)*Dxx) + ...
    1i/2*(spdiags(y3pp, 0, my3, my3)*Dx); % New operator for 3rd dimension

%% Vortice cilindrico infinito (entro i limiti della griglia scelti dall'utente)
X1 = antitransform(y1); 
X2 = antitransform(y2);
X3 = antitransform(y3);

U0 = sqrt(R).*exp(1i*atan2(X2,X1));
%% Strang time integration 
U = U0;

tau = 0.05;
tstar = 20;
t = 0;

E1 = expm(tau*A1);
E2 = expm(tau*A2);
E3 = expm(tau*A3); 

while t < tstar
    
    U = ndcovlt(U,{E1,E2,E3});
    U = exp(tau*(1i/2)*(1 - abs(U).^2)).*U; 
    U = ndcovlt(U,{E1,E2,E3});
    
    t = t + tau;
    disp(['Time: ', num2str(t), ' / ', num2str(tstar)]); % Display progress
end

%% plot final solution (adapt for 3D visualization)
% Visualizing 3D data can be challenging. We can plot slices or isosurfaces.

figure
val = 0.5; % livello di densità
p = patch(isosurface(X,Y,Z,abs(U).^2,val));
isonormals(X,Y,Z,abs(U).^2,p);
set(p,'FaceColor','red','EdgeColor','none','FaceAlpha',0.5);
daspect([1 1 1]); camlight; lighting gouraud;
xlabel('x'); ylabel('y'); zlabel('z');
title('|U|^2 finale (isosurface)');

