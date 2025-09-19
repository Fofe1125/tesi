function m = cmass(rho, X, Y)

% integral over domain X,Y of rho
% using Chebyshev-Gauss quadrature
% [X,Y] = ndgrid(x,y);

m = sum(rho(:)).*(X(1,2) - X(1,1)).*(Y(2,1) - Y(1,1));

end
