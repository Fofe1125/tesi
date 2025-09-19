function R = rho3d(x, y, z, curve_points)
% rho3D_curve - 4-th order Padè approximation of a 3D vortex around a parameterized curve.
%
%   R = rho3d(x, y, z, curve_points) calculates the density R at points
%   (x, y, z) based on their minimum distance to a 3D parameterized curve.
%
%   Inputs:
%     x, y, z        - Arrays (same size) representing the coordinates where
%                      the density R is to be evaluated.
%     curve_points   - A Nx3 matrix where each row [xp, yp, zp] represents
%                      a point on the parameterized curve. The curve is
%                      approximated by these discrete points.
%
%   Output:
%     R              - Array (same size as x,y,z) with the density values.

% Padè approximation coefficients (as provided for the 2D case)
a(1) = 0.34010790700196714760;
b(1) = (2304*a(1)^3 +656*a(1)^2 -421*a(1) - 28)/(7680*a(1)^2 - 1680*a(1) - 330);
a(2) = a(1)*(b(1) - 1/4);
b(2) = ((737280*a(1)^3 + 209920*a(1)^2 - 134720*a(1) - 8960)*b(1) - ...
    364544*a(1)^3 + 70144*a(1)^2 + 18256*a(1) + 393)/( ...
    2457600*a(1)^2 - 537600*a(1) - 105600);
a(3) = (a(1)*(192*b(2) - 48*b(1) + 16*a(1) + 5))/192;
b(3) = ((61440*a(1) - 9600)*b(2) + (-30720*a(1)^2 + 640*a(1) + 560)*b(1) + ...
    8448*a(1)^2 - 1056*a(1) - 21)/(368640*a(1) - 92160);
a(4) = (4608*a(1)*b(3) - 1152*a(1)*b(2) + (384*a(1)^2 + 120*a(1))*b(1) - ...
    128*a(1)^2 - 7*a(1))/4608;

% preallocate
R = zeros(size(x));

xv = x(:);
yv = y(:);
zv = z(:);

num_eval_points = numel(xv);
min_dist_sq = zeros(num_eval_points, 1);

for i = 1:num_eval_points
    dx = xv(i) - curve_points(:,1);
    dy = yv(i) - curve_points(:,2);
    dz = zv(i) - curve_points(:,3);
    distances_to_curve_sq = dx.^2 + dy.^2 + dz.^2;
    min_dist_sq(i) = min(distances_to_curve_sq);
end

Rvals = (a(1).*min_dist_sq + a(2).*min_dist_sq.^2 + ...
         a(3).*min_dist_sq.^3 + a(4).*min_dist_sq.^4) ./ ...
        (1 + b(1).*min_dist_sq + b(2).*min_dist_sq.^2 + ...
         b(3).*min_dist_sq.^3 + a(4).*min_dist_sq.^4);

R = reshape(Rvals, size(x));

end
