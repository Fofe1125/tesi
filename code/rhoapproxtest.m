% solve_compactified.m
clearvars; close all; clc

% discretizzazione su s in [0,1]
m = 501;
s = linspace(0,1,m)';      % s(1)=0, s(end)=1
h = s(2)-s(1);

% derivata prima e seconda (matrici centrali)
e = ones(m,1);
D1 = spdiags(e*[-1,0,1]/(2*h), -1:1, m,m);
D2 = spdiags(e*[1,-2,1]/(h^2), -1:1, m,m);

% gestione bordi: per ora non usiamo le righe 1 e m delle matrici (imporremo BC)
D1(1,:) = 0; D1(m,:) = 0;
D2(1,:) = 0; D2(m,:) = 0;

% funzione residuo f(g) (vector), e Jacobiano Jf(g)
% NOTA: usare il segno corretto (1 - g.^2).*g e Jacobiano 1 - 3 g.^2
nonlinear_term = @(g) (1 - g.^2).*g;
nonlinear_Jdiag = @(g) (1 - 3*g.^2);   % derivata di (1-g^2)g rispetto a g

Ffun = @(g) ( (s-1).^4 ).* (D2*g) + 2*(s-1).^3 .* (D1*g) - (((s-1).^3)./s).* (D1*g) ...
           - (((s-1).^2)./(s.^2)).* g + nonlinear_term(g);

Jfun = @(g) spdiags((s-1).^4,0,m,m)*D2 + spdiags(2*(s-1).^3,0,m,m)*D1 ...
           - spdiags(((s-1).^3)./s,0,m,m)*D1 - spdiags(((s-1).^2)./(s.^2),0,m,m) ...
           + spdiags(nonlinear_Jdiag(g),0,m,m);

% condizioni al contorno: g(1)=0, g(end)=1
% costruiamo il sistema su tutti gli m punti ma sovrascriviamo le equazioni di bordo
% con condizioni essenziali (Dirichlet).
g = (linspace(0,1,m)').^2;    % iniziale
g(1) = 0;
g(end) = 1;

tol = 1e-8;
maxit = 100;
for k = 1:maxit
    F = Ffun(g);
    J = Jfun(g);
    % impongo le condizioni di Dirichlet:
    F(1) = g(1) - 0;
    F(end) = g(end) - 1;
    % riga identità nel Jacobiano per le BC
    J(1,:) = 0; J(1,1) = 1;
    J(end,:) = 0; J(end,end) = 1;
    
    % risolvo Newton linear system
    delta = - (J \ F);
    
    g = g + delta;
    g(1) = 0; g(end) = 1;   % reinforzo BC (numerical safety)
    
    if norm(delta, Inf) < tol
        fprintf('Newton converged in %d iters, ||delta||_inf = %.3e\n', k, norm(delta,Inf));
        break;
    end
    if k==maxit
        warning('Newton did not converge within maxit');
    end
end

% Ricostruzione r, f e due candidate per rho
s_interior = s;                 % 1..m (include 0 and 1)
r = s_interior ./ (1 - s_interior);   % r(1)=0, r(end)=Inf (NaN at s=1)
r(end) = 1e6; % evita Inf nella vettorizzazione (sostituisci con un grande numero)
f = g;                          % per definizione g(s)=f(r)

% DUE possibili formule per rho (scegli quella corretta per la tua definizione)
% Variante A (esempio): rho_A = ( g ./ (1 + r) ).^2;
rho_A = ( g ./ (1 + r) ).^2;

% Variante B (altro esempio): rho_B = ( (r .* g) ./ (1 + r) ).^2;
rho_B = ( (r .* g) ./ (1 + r) ).^2;

% Correggi valore in r=0 (i=1): usare limite se necessario (qui usiamo il valore ottenuto)
rho_A(1) = 0;
rho_B(1) = 0;

% visualizzo risultati
figure;
subplot(2,2,1), plot(r, g,'.-'); xlabel('r'); title('f(r) = g(s)');
subplot(2,2,2), plot(r, rho_A,'.-'); xlabel('r'); title('\rho_A(r) = (g/(1+r))^2');
subplot(2,2,3), plot(r, rho_B,'.-'); xlabel('r'); title('\rho_B(r) = (r g/(1+r))^2');
subplot(2,2,4), plot(s, g,'.-'); xlabel('s'); title('g(s) on compact domain');

% salva risultati (r, f, rho)
res.r = r;
res.f = f;
res.rho_A = rho_A;
res.rho_B = rho_B;
save('solution_compactified.mat','res','s','g');

fprintf('Salvati risultati in solution_compactified.mat (variabili: res, s, g)\n');
