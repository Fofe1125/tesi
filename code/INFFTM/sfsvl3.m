function psipdb = sfsvl3(point,vector,x)
% Function file: psipdb = sfsvl3(point,vector,x)
%
% Description: Compute psi originating from a single straight 3d vortex.
%
% Input variables:
%   point:  a three-dimensional point, as a row vector, in the core of a vortex
%   vector: an oriented direction of a vortex, as a row vector
%   x:      cell array containing the column vectors x{d}, d=1:3, of the
%           coordinates of the points in the physical domain. From x
%           it is possible to generate the grid points [X{1:3}] = ndgrid(x{1:3})
%
% Output variables:
%   psipdb: complex 3d-array of the solution in the physical domain
[X{1:3}] = ndgrid(x{1:3});
eta = size(X{1});
for d = 1:3
  X{d} = X{d}(:);
end
% choose a meaningful normal vector (nv) and normalize the tangent
% vector (tv)
tv = vector/norm(vector);
nv = zeros(size(tv));
[~,idx] = sort(abs(tv));
nv(idx(1)) = sign(tv(idx(3)));
nv(idx(2)) = 0;
nv(idx(3)) = -tv(idx(1))/tv(idx(3));
nv = nv/norm(nv);
% compute the binormal
bv = cross(tv,nv);
% For each grid point, tbar is the value of the parameter
% corresponding to minimal distance from the line, that is
% (x-point-tbar*tv)*tv=0
tbar = ([X{1},X{2},X{3}]-repmat(point,length(X{1}),1))*tv(:);
% vectors from point+tbar*tv to the corresponding grid point
rvec = [X{1},X{2},X{3}]-...
       (repmat(point,length(X{1}),1)+repmat(tbar,1,3)*diag(tv));
% Components of the vector rvec along the binormal and normal
DB = rvec*bv(:); DN = rvec*nv(:);
% Compute phase
theta = reshape(atan2(DB(:),DN(:)),eta);
% Compute density according to Pade' approximation
r2 = reshape(DB.^2+DN.^2,eta);
psipdb = sqrt(sf4pade(r2)).*exp(1i*theta);
