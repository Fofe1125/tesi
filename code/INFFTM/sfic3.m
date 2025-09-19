function sf = sfic3(x,pdb,mirror,point,vector,t0)
% Function file: sf = sfic3(x,pdb,mirror,point,vector)
% Function file: sf = sfic3(x,pdb,mirror,point,vector,t0)
%
% Description: Compute the structure sf for superimposed 3d straight vortices.
%
% Input variables:
%   x:      cell array containing the column vectors x{d}, d=1:3, of the
%           coordinates of the points in the physical domain. From x it is
%           possible to generate the grid points [X{1:3}] = ndgrid(x{1:3})
%   pdb:    physical domain boundaries in the axis form
%   mirror: row vector of mirroring flags true/false
%   point:  2d-array identifying, for each row, a three-dimensional
%           point in the core of a vortex
%   vector: 2d-array identifying, for each row, an oriented direction
%           of a vortex
%  t0:      optional initial time (default is 0).
%
%  Output variables:
%   sf:     structure defining the function psi, where
%           sf.pdb:    are the physical domain boundaries
%           sf.psipdb: is the complex 3d-array of psi in the physical domain
%           sf.mirror: is the row vector of mirroring flags true/false
%           sf.t:      is the initial time (usually, but not necessarily, 0)
sf.psipdb = sfsvl3(point(1,:),vector(1,:),x);
for i = 2:size(point,1)
  sf.psipdb = sf.psipdb .* sfsvl3(point(i,:),vector(i,:),x);
end
sf.pdb = pdb;
if (nargin == 5)
  t0 = 0;
end
sf.t = t0;
sf.mirror = mirror;
