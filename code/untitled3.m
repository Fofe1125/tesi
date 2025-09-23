clear all;
close all;

L = 30;
mx = 347;
my = 347;
hx = 2*L/(mx - 1);
hy = 2*L/(my - 1);

Dxx = spdiags(ones(mx,1)*[1, -2, 1]/(hx^2), -1:1, mx,mx);
Dyy = spdiags(ones(my,1)*[1, -2, 1]/(hy^2), -1:1, my,my);

% Homogeneous Neumann
Dxx(1,1:2) = [-2,2]/(hx^2);
Dyy(1,1:2) = [-2,2]/(hy^2);
Dxx(mx,mx-1:mx) = [2,-2]/(hx^2);
Dyy(my,my-1:my) = [2,-2]/(hy^2);

x = linspace(-L, L, mx)';
y = linspace(-L, L, my)';

[X, Y] = ndgrid(x, y);

%% initial condition
rhoVortex = rho(X, Y, 0,0); % definisci la tua funzione rho
phase = atan2(Y,X);

U0 = sqrt(rhoVortex).*exp(1i*phase);
u0 = reshape(U0,mx*my,1);

%% Implicit time integration (Newton-Krylov)
u = u0;
t = 0;
k = 0.1; % time step
maxTime = 1;
tolNewton = 1e-6;

while t < maxTime
    u_old = u;
    
    % Definizione funzione non lineare F(u) = 0
    F = @(u) u - u_old - 1i/2*k*pmi_step((u + u_old)/2, Dxx, Dyy, mx, my);

    % Newton-Krylov iterations
    maxNewton = 20;
    for newtonIter = 1:maxNewton
        res = F(u);
        normRes = norm(res);
        fprintf('t = %.3f, Newton iter %d, ||F|| = %.3e\n', t, newtonIter, normRes);

        if normRes < tolNewton
            break;
        end

        % funzione anonima per J*v
        Jv = @(v) jacobian_vector_product(F, u, v);

        % Risolvi J * delta = -F(u) usando GMRES
        [delta, flag] = gmres(Jv, -res, [], 1e-6, 50); % restart 50
        if flag ~= 0
            warning('GMRES did not converge');
        end

        u = u + delta;
    end

    t = t + k;
end

U = reshape(u, mx, my);

function v = pmi_step(x, Dxx, Dyy, mx, my)
% PMI_STEP  calcola l'operatore PMI step su x (vettore mx*my)
%
% x   : vettore colonna (mx*my)
% Dxx : derivata seconda lungo x (mx x mx)
% Dyy : derivata seconda lungo y (my x my)
% mx, my : dimensioni griglia

    % Riscalare x in matrice 2D
    Xmat = reshape(x, mx, my);

    % Applicare derivata secondo x e y usando prodotti Kronecker
    % Dxx * X + X * Dyy' equivale a (Dxx kron I + I kron Dyy) * vec(X)
    Lx = Dxx * Xmat;      % mx x my
    Ly = Xmat * Dyy.';    % mx x my

    % Termine non lineare locale: (1 - |x|^2)*x
    Nl = (1 - abs(Xmat).^2) .* Xmat;

    % Ricombina tutto e rimette in vettore colonna
    Vmat = Lx + Ly + Nl;
    v = reshape(Vmat, mx*my, 1);
end

