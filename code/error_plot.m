clearvars;
close all;

errN303011001 = csvread("data/RE_NUG_L30_301S_1000T_FD.csv");
errU303011001 = csvread("data/RE_UG_L30_301S_1000T_FD.csv");
errFB = csvread("data/error_tfd4.csv");

tspan = errN303011001(:,1);

errN = errN303011001(:,2);
errU = errU303011001(:,2);
err4 = errFB(:,2);

figure;
hold on;
set(gca,'YScale','log');
plot(tspan, errU, '-rx');
plot(tspan, errN, '-g^');
plot(tspan, err4, '-bd');
xlabel('Time (s)');
ylabel('Error Magnitude');
legend('UFD, L = 30, m = 301','FD, L = 30, m = 301', 'FBTD, m = 87');
title('Error Comparison');

% etichette assi
xlabel('t', 'FontSize', 12);
ylabel('Max pointwise relative error', 'FontSize', 12);
