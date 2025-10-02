clear ,all;
close all;

load('data/grid_ufd.mat', 'X', 'Y');
X_ufd = X; Y_ufd = Y;

load('data/grid_nfd.mat', 'X', 'Y');
X_nfd = X; Y_nfd = Y;

%% grid plot
% Definisci i nodi della griglia

L = 5;

figure;
subplot(1,2,1);
hold on;
% Linee verticali (costanti in x)
plot(X_ufd', Y_ufd', 'k-', 'LineWidth', 0.1);

% Linee orizzontali (costanti in y)
plot(X_ufd, Y_ufd, 'k-', 'LineWidth', 0.1);

xlabel('x'); ylabel('y');
title('Uniform grid');
xlim([-L, L]);
ylim([-L, L]);
hold off;

subplot(1,2,2);
hold on;
% Linee verticali (costanti in x)
plot(X_nfd', Y_nfd', 'k-', 'LineWidth', 0.1);

% Linee orizzontali (costanti in y)
plot(X_nfd, Y_nfd, 'k-', 'LineWidth', 0.1);

xlabel('x'); ylabel('y');
title('Non uniform grid');
xlim([-L, L]);
ylim([-L, L]);
hold off;