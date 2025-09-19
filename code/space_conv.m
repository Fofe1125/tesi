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
m = 501;
h = 2/(m - 1);
y1 = linspace(-1 + epsbd,1 - epsbd,m);
y2 = linspace(-1 + epsbd,1 - epsbd,m);

y1p = 2./(pi*alpha*(tan(pi/2*y1).^2 + 1));
y2p = 2./(pi*alpha*(tan(pi/2*y2).^2 + 1));
y1pp = -(4*tan(pi/2*y1))./(pi*alpha^2*(tan(pi/2*y1).^2 + 1).^2);
y2pp = -(4*tan(pi/2*y2))./(pi*alpha^2*(tan(pi/2*y2).^2 + 1).^2);

[Y1,Y2] = ndgrid(y1,y2);

% laplacian in mapped space
[Dx, Dxx] = buildMatrix(m, h);
A1 = 1i/2*(spdiags(y1p.^2, 0, m, m))*Dxx + ...
    1i/2*(spdiags(y1pp, 0, m, m))*Dx;
A2 = 1i/2*(spdiags(y2p.^2, 0, m, m))*Dxx + ...
    1i/2*(spdiags(y2pp, 0, m, m))*Dx;

%% vortex 
X1 = antitransform(Y1); X2 = antitransform(Y2);
% rhoVortex = rho(Y1, Y2, 0, 0);
% phaseVortex = atan2(Y2, Y1);
rhoVortex = rho(X1, X2, 0,0);
phaseVortex = atan2(X2,X1);

U0 = sqrt(rhoVortex).*exp(1i*phaseVortex);

%% time splitting integration
tstar = 0.1;
U = U0;
nt = 100;
tau = tstar/(nt - 1);

E1 = expm(tau/2*A1);
E2 = expm(tau/2*A2);

for k = 1:nt - 1
    
    U = E1*U*E2.';
    U = exp(tau*(1i/2)*(1 - abs(U).^2)).*U;
    U = E1*U*E2.';
 
end

Uref = U;
mref = m;

mrange = 10:10:100;
counter = 1;
err = NaN(1, length(mrange));

for m = mrange + 1

    fprintf("Starting sim with m = %d\n", m);

    h = 2/(m - 1);
    y1 = linspace(-1 + epsbd,1 - epsbd,m);
    y2 = linspace(-1 + epsbd,1 - epsbd,m);
    
    y1p = 2./(pi*alpha*(tan(pi/2*y1).^2 + 1));
    y2p = 2./(pi*alpha*(tan(pi/2*y2).^2 + 1));
    y1pp = -(4*tan(pi/2*y1))./(pi*alpha^2*(tan(pi/2*y1).^2 + 1).^2);
    y2pp = -(4*tan(pi/2*y2))./(pi*alpha^2*(tan(pi/2*y2).^2 + 1).^2);
    
    [Y1,Y2] = ndgrid(y1,y2);
    
    % laplacian in mapped space
    [Dx, Dxx] = buildMatrix(m, h);
    A1 = 1i/2*(spdiags(y1p.^2, 0, m, m))*Dxx + ...
        1i/2*(spdiags(y1pp, 0, m, m))*Dx;
    A2 = 1i/2*(spdiags(y2p.^2, 0, m, m))*Dxx + ...
        1i/2*(spdiags(y2pp, 0, m, m))*Dx;
    
    X1 = antitransform(Y1); X2 = antitransform(Y2);
    rhoVortex = rho(X1, X2, 0,0);
    phaseVortex = atan2(X2,X1);
    
    U0 = sqrt(rhoVortex).*exp(1i*phaseVortex);
    U = U0;
    
    E1 = expm(tau/2*A1);
    E2 = expm(tau/2*A2);
    
    for k = 1:nt - 1
        
        U = E1*U*E2.';
        U = exp(tau*(1i/2)*(1 - abs(U).^2)).*U;
        U = E1*U*E2.';
     
    end

    mask = (abs(U0) > 0);

    N = abs(U - Uref(1:(mref - 1)/(m - 1):mref, 1:(mref - 1)/(m - 1):mref));
    D =  abs(Uref(1:(mref - 1)/(m - 1):mref, 1:(mref - 1)/(m - 1):mref));
    relative_error = N(mask) ./ D(mask);
    err(counter) = norm(relative_error(:), 'inf');
    counter = counter + 1;

end

%% space conv
figure;
set(gca,'Xscale', 'log', 'Yscale', 'log');
hold on;
plot(mrange, err, '--bd');
plot(mrange,err(end)*(mrange/mrange(end)).^(-2),'m');
xlabel('Time steps');
ylabel('Pointwise error')
legend('Pointwise error', '2nd order', 'Location', 'southwest');
grid on;