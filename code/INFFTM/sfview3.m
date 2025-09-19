function x = sfview3(sf,iso)
% Function file: x = sfview3(sf,iso)
%
% Description: Show an isosurface of the density associated to a structure sf.
%
% Input variables:
%   sf:     structure defining the function psi, where
%           sf.pdb:    are the physical domain boundaries
%           sf.psipdb: is the complex 3d-array of psi in the physical domain
%           sf.mirror: is the row vector of mirroring flags true/false
%           sf.t:      is the simulation time
%   iso:    isosurface level of the density
%
% Output variables:
%   x:      cell array containing the column vectors x{d} of the
%           coordinates of the points in the physical domain.
%           From x it is possible to generate the grid points
%           [X{1:end}] = ndgrid(x{1:end})
eta = size(sf.psipdb);
x = cell(1,3);
for d = 1:3
  x{d} = linspace(sf.pdb(2*d-1),sf.pdb(2*d),eta(d)+1)';
  x{d} = x{d}(1:eta(d));
end
plotiso3(x,real(sf.psipdb).^2+imag(sf.psipdb).^2,iso)
title(sprintf('simulation time = %g, isosurface rho = %g',sf.t,iso))
