function plotiso3(x,data,iso)
% Function file: plotiso3(x,data,iso)
%
% Description: Plot an isosurface of the real 3d input data.
%
% Input variables:
%   x:      cell array containing the column vectors x{d}, d=1:3, of the
%           coordinates of the points in the physical domain. From x it is
%           possible to generate the grid points [X{1:3}] = ndgrid(x{1:3})
%   data:   real data corresponding to ndgrid(x{1:3})
%   iso:    isosurface level to display

% Matlab requires the input of isosurface to be in meshgrid format
data = permute(data,[2,1,3]);
[X{1:3}] = meshgrid(x{1:3});
figure
isosurface(X{1},X{2},X{3},data,iso);
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
box on
