clearvars;
close all;

%% load data
load('err_nfd.mat', 'times', 'err');
times_nfd = times; err_nfd = err;

load('err_ufd.mat', 'times', 'err');
times_ufd = times; err_ufd = err;

%% plot results
figure;
hold on;
set(gca,'Yscale', 'log');
plot(times_nfd, err_nfd, 'bd--', 'DisplayName', 'NFD Error');
plot(times_ufd, err_ufd, 'r^--', 'DisplayName', 'UFD Error');
xlabel('Time');
ylabel('Error');
title('Error Comparison');
legend show;
grid on;