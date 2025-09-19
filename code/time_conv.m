% Neumann Bc truncated domain finite diff 2-nd order
% stationary vortex

clearvars;
close all;
%% space discretization
tic;
L = 30;
Nx = 301;
Ny = 301;
hx = 2*L/(Nx - 1);
hy = 2*L/(Ny - 1);

Dxx = spdiags(ones(Nx,1)*[1, -2, 1]/(hx^2), -1:1, Nx,Nx);
Dyy = spdiags(ones(Ny,1)*[1, -2, 1]/(hy^2), -1:1, Ny,Ny);

% Homogeneus Neumann conditions 
Dxx(1,1:2) = [-2,2]/(hx^2);
Dyy(1,1:2) = [-2,2]/(hy^2);
Dxx(Nx,Nx-1:Nx) = [2,-2]/(hx^2);
Dyy(Ny,Ny-1:Ny) = [2,-2]/(hy^2);

epsbd = 1e-8;
x = linspace(-L+epsbd, L-epsbd, Nx)';
y = linspace(-L+epsbd, L-epsbd, Ny)';

[X, Y] = ndgrid(x, y);

%% initial condition
rhoVortex = rho(X, Y, 0,0);
phase = atan2(Y,X);

U0 = sqrt(rhoVortex).*exp(1i*phase);
% error measure
Rgrid = sqrt(X.^2 + Y.^2);
R = 100;       
tol = 1e-6;           
mask = (Rgrid < R & abs(U0) > tol); % avoid division by 0
%% time splitting integration - ref solution
tstar = 0.1;
ntref = 1500;

tau = tstar/ntref;

Ex = expm(1i*(tau/4)*Dxx);
Ey = expm(1i*(tau/4)*Dyy);

U = U0;

for k = 1:ntref 

    U = Ex*U*Ey.';
    U = exp(tau*(1i/2)*(1 - abs(U).^2)).*U;
    U = Ex*U*Ey.';

end

Uref = U;    

%% temporal convergence
ntrange = 10:10:300;
taurange = tstar./ntrange;
relative_error = NaN(1,length(taurange));

parfor i = 1:length(taurange)

    U = U0;
    tau = taurange(i);
    fprintf("Starting simulation with tau = %.5f\n", tau);

    Ex = expm(1i*tau/4*Dxx);
    Ey = expm(1i*tau/4*Dyy);

    for k = 1:ntrange(i) 
        % time step
        U = Ex*U*Ey.';
        U = exp(tau*(1i/2)*(1 - abs(U).^2)).*U;
        U = Ex*U*Ey.';
         
    end

    % error measure
    re_err = abs(U(mask) - Uref(mask)) ./ abs(Uref(mask));
    relative_error(i) = norm(re_err(:), 'inf');
    fprintf('Relative error at tau = %.5f:%.16f\n', tau, relative_error(i));
end
toc;
%% convergence order 
figure;
set(gca,'Xscale', 'log', 'Yscale', 'log');
hold on;
plot(ntrange, relative_error, '--bd');
plot(ntrange,relative_error(end)*(ntrange/ntrange(end)).^(-2),'m');
xlabel('Time steps');
ylabel('Pointwise error')
legend('Pointwise error', '2nd order', 'Location', 'southwest');
grid on;

