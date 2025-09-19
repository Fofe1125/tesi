% debug_vortex_mapped.m  -- versione diagnostica da copiare/incollare
clearvars; close all; clc;

%% space rescalation
alpha = 5;
transform = @(xl) (2/pi)*atan(xl/alpha);
antitransform = @(yl) alpha*tan((pi/2)*yl);

%% space discretization of (-1,1)^2
epsbd = 1e-8;
my1 = 201; my2 = 201;
h = 2/(my1 - 1);
y1 = linspace(-1 + epsbd,1 - epsbd,my1);
y2 = linspace(-1 + epsbd,1 - epsbd,my2);

% derivatives of map
y1p = 2./(pi*alpha*(tan(pi/2*y1).^2 + 1));
y2p = 2./(pi*alpha*(tan(pi/2*y2).^2 + 1));
y1pp = -(4*tan(pi/2*y1))./(pi*alpha^2*(tan(pi/2*y1).^2 + 1).^2);
y2pp = -(4*tan(pi/2*y2))./(pi*alpha^2*(tan(pi/2*y2).^2 + 1).^2);

[Y1,Y2] = ndgrid(y1,y2);    % attenzione all'orientamento: Y1 varia su righe

% laplacian in mapped space
[Dx, Dxx] = buildMatrix(my1, h); % verifica che buildMatrix dia gli scaling corretti
A1 = 1i/2*(spdiags(y1p.^2, 0, my1, my1))*Dxx + 1i/2*(spdiags(y1pp, 0, my1, my1))*Dx;
A2 = 1i/2*(spdiags(y2p.^2, 0, my2, my2))*Dxx + 1i/2*(spdiags(y2pp, 0, my2, my2))*Dx;

%% vortex nel dominio fisico, valutato sui nodi mappati
X1 = antitransform(Y1); X2 = antitransform(Y2);

xc = 3; yc = 0;
rhoVortex = rho(X1, X2, xc,yc).*rho(X1,X2,-xc,yc);
phaseVortex = atan2(X2 - yc, X1 - xc) - atan2(X2 + yc, X1 + xc);
U0 = sqrt(rhoVortex).*exp(1i*phaseVortex);

%% controlli iniziali rapidi
fprintf('X1 range: [%g, %g], X2 range: [%g, %g]\n', min(X1(:)), max(X1(:)), min(X2(:)), max(X2(:)));
fprintf('|U0| min/max = %g / %g\n', min(abs(U0(:))), max(abs(U0(:))));

% Jacobiano per integrali fisici
J = (y1p(:))*(y2p(:)).';
mass0_phys = h*h * sum(abs(U0(:)).^2 .* J(:));
fprintf('Mass (phys) initial = %.6e\n', mass0_phys);

%% diagnostica operatori
% norma stimata degli operatori lineari
try
    nA1 = normest(A1); nA2 = normest(A2);
catch
    nA1 = norm(full(A1)); nA2 = norm(full(A2));
end
fprintf('Estimated norm A1 = %.3e, A2 = %.3e\n', nA1, nA2);

%% test propagatori con tau piccolo
tstar = 100;
nt = 200;               % prova con passo moderato
tau = tstar/(nt - 1);

E1 = expm(tau/2*A1);
E2 = expm(tau/2*A2);

% quanto sono vicini all'identità?
I1 = eye(size(E1)); I2 = eye(size(E2));
fprintf('max |E1-I| = %.3e, max |E2-I| = %.3e\n', max(abs(E1(:)-I1(:))), max(abs(E2(:)-I2(:))));

% applica singolo passo lineare test
Utest = E1*U0*E2.';
fprintf('max |Utest - U0| = %.3e\n', max(abs(Utest(:) - U0(:))));

% se la differenza è ~0: o A1/A2 ~ 0 (scaling errato) o tau troppo piccolo
if max(abs(Utest(:)-U0(:))) < 1e-12
    warning('Lineare applicazione produce praticamente identita: controlla scaling Dxx/Dx o tau.');
end

%% funzione centro di massa (con jacobiano)
compute_com = @(U) deal( ...
    sum( (abs(U).^2 .* J) .* X1 , 'all') / sum(abs(U).^2 .* J,'all'), ...
    sum( (abs(U).^2 .* J) .* X2 , 'all') / sum(abs(U).^2 .* J,'all') ...
    );

[cx0, cy0] = compute_com(U0);
fprintf('Center of mass initial: (%.6f, %.6f)\n', cx0, cy0);

%% prova: piccola perturbazione (rompe eventuale autostazionarietà)
rng(1);
pert = 1 + 1e-3*(randn(size(U0)) + 1i*randn(size(U0)));
U = U0 .* pert;

%% integrazione in tempo con diagnostica
nt = 2000;            % aumenta nt per sicurezza
tstar = 100;
tau = tstar/(nt - 1);
E1 = expm(tau/2*A1); E2 = expm(tau/2*A2);

COM = zeros(nt,2);
COM(1,:) = [cx0, cy0];
mass_track = zeros(nt,1);
mass_track(1) = h*h * sum(abs(U(:)).^2 .* J(:));

fprintf('Starting time loop nt=%d, tau=%.3e\n', nt, tau);
for k = 1:nt-1
    Uprev = U;
    U = E1*U*E2.';                                     % passo lineare
    U = exp(tau*(1i/2)*(1 - abs(U).^2)).*U;             % passo nonlineare (fase)
    U = E1*U*E2.';                                     % passo lineare
    
    % diagnosi ogni 10%
    if mod(k,round(nt/10))==0 || k==1
        dmax = max(abs(U(:)-Uprev(:)));
        [cx,cy] = compute_com(U);
        fprintf('step %d/%d: max|dU|=%.3e, COM=(%.6f,%.6f)\n', k, nt-1, dmax, cx, cy);
    end
    
    COM(k+1,:) = compute_com(U);
    mass_track(k+1) = h*h * sum(abs(U(:)).^2 .* J(:));
end

fprintf('Mass relative change = %.3e\n', (mass_track(end)-mass_track(1))/mass_track(1));

%% plot risultati finali e traiettoria COM
lim = 15; L = [-lim lim];
figure('Name','Results');
subplot(2,2,1);
surf(X1,X2,abs(U0),'EdgeColor','none'); view(2); title('|U_0|'); colorbar; xlim(L); ylim(L);

subplot(2,2,2);
surf(X1,X2,abs(U),'EdgeColor','none'); view(2); title('|U final|'); colorbar; xlim(L); ylim(L);

subplot(2,2,3);
surf(X1,X2,angle(U),'EdgeColor','none'); view(2); title('phase final'); colorbar; xlim(L); ylim(L);

subplot(2,2,4);
plot(linspace(0,tstar,nt), COM(:,1), 'LineWidth',1.2); hold on;
plot(linspace(0,tstar,nt), COM(:,2), 'LineWidth',1.2); legend('COM_x','COM_y'); xlabel('t'); title('COM trajectory');

% mostra differenza assoluta
figure('Name','difference abs');
surf(X1,X2,abs(U)-abs(U0),'EdgeColor','none'); view(2); title('abs(U)-abs(U0)'); colorbar; xlim(L); ylim(L);
