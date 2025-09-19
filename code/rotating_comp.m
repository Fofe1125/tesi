clear all;
close all;

% comparison between ufd and nfd on rotating vortex
load('data/rotating_nfd.mat', 'U', 'X', 'Y');
U_nfd = U; X_nfd = X; Y_nfd = Y;

load('data/rotating_ufd.mat', 'U', 'X', 'Y');
U_ufd = U; X_ufd = X; Y_ufd = Y;

%% plot results
figure;

subplot(1,2,1);
surf(X_ufd, Y_ufd, abs(U_ufd), 'EdgeColor','none');
colorbar; colormap('jet');
axis equal tight;
title('UFD, m = 771')
view(2);

subplot(1,2,2);
surf(X_nfd, Y_nfd, abs(U_nfd), 'EdgeColor','none');
colorbar; colormap('jet');
axis equal tight;
title('NFD, m = 255')
view(2);

