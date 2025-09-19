function x = mkgrid(pdb,eta)
% Function file: x = mkgrid(pdb,eta)
%
% Description: Generate a rectilinear grid in the physical domain.
%
% Last point(s) in the physical domain are removed by periodicity.
%
% Input variables:
%   pdb:    physical domain boundaries in the axis form
%   eta:    row vector of numbers of grid points in the physical domain,
%           eta(d) even
%
% Output variables:
%   x:      cell array containing the column vectors x{d} of the
%           coordinates of the points in the physical domain.
%           From x it is possible to generate the grid points
%           [X{1:end}] = ndgrid(x{1:end})
dd = length(eta);
x = cell(1,dd);
for d = 1:dd
  x{d} = linspace(pdb(2*d-1),pdb(2*d),eta(d)+1)';
  x{d} = x{d}(1:eta(d)); % remove last point by periodicity
end
%!test
%! x = mkgrid([0,1],10);
%! assert(x{1},linspace(0,1,11)(1:10)')
%!test
%! x = mkgrid([0,1,-1,1],[10,10]);
%! assert(x{2},linspace(-1,1,11)(1:10)')
%!test
%! x = mkgrid([0,1,-1,1,-10,10],[10,10,20]);
%! assert(x{3},linspace(-10,10,21)(1:20)')
%!demo
%! x = mkgrid([0,1,-1,1,-10,10],[10,10,20]);
