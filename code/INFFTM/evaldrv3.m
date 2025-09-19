% Description: Script to perform fast evaluation at grids or artibrary points.

% Before running this script, make sure NFFT is available and set the
% correct path in nfftpath.m.
axislim = [-5,5,-4,6,-3,7];
load 'data00.mat'
iso = 0.1;
sfview3(sf,iso);
axis(sf.pdb)
view(8,16)
drawnow
%print('-depsc','rho0.eps')
load 'data10.mat'
iso = 0.1;
sfview3(sf,iso);
% plot of orignal data
axis(sf.pdb)
view(8,16)
drawnow
%print('-depsc','rho10.eps')
iso = 0.05;
sfview3(sf,iso);
view(8,16)
axis(axislim)
drawnow
%print('-depsc','rho10zoom.eps')
% evaluation at the rectilinear equispaced grid,
% number of evaluation points in [a(1),b(1)]x[a(2),b(2)],[a(3),b(3)]
eta = size(sf.psipdb);
M = 8*eta;
y = cell(1,3);
for d = 1:3
  y{d} = linspace(sf.pdb(2*d-1),sf.pdb(2*d),M(d)+1)';
end
[psihat,a,b] = sf2psihat(sf);
psi = igridftn(psihat,a,b,y);
iso = 0.0012;
plotiso3(y,real(psi).^2+imag(psi).^2,iso)
title(sprintf('simulation time = %g, isosurface rho = %g',sf.t,iso))
view(8,16)
axis(axislim)
drawnow
%print('-depsc','rho10eq.eps')
% evaluation at the rectilinear non-equispaced grid
M = 7*eta;
k = [1.02,1.02,1.02];
for d = 1:3
  h0 = sf.pdb(2*d)*(k(d)-1)/(k(d)^(M(d)/2)-1);
  y{d} = h0*cumsum(k(d).^(0:M(d)/2-1)');
  y{d} = [-flipud(y{d});0;y{d}];
end
psi = igridftn(psihat,a,b,y);
iso = 0.0012;
plotiso3(y,real(psi).^2+imag(psi).^2,iso)
title(sprintf('simulation time = %g, isosurface rho = %g',sf.t,iso))
view(8,16)
axis(axislim)
drawnow
%print('-depsc','rho10neq.eps')
% nfft-evaluation
[Xi,rho] = sftubeeval3(sf,[0.2,0.05]);
iso = 0.0012;
ind = find(abs(rho <= iso));
figure
plot3(Xi(1,ind),Xi(2,ind),Xi(3,ind),'.k','MarkerSize',4)
box on
title(sprintf('simulation time = %g, isosurface rho = %g',sf.t,iso))
axis equal
axis(axislim)
view(8,16)
drawnow
xlabel('x')
ylabel('y')
zlabel('z')
%print('-depsc','rho10nfft.eps')
