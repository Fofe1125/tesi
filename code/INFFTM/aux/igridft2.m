function psi = igridft2(psihat,a,b,y)
% Function file: psi = igridft2(psihat,a,b,y)
%
% Description: Evaluate 2-dimensional Fourier coeffs to a rectilinear grid.
%
% Input variables:
%   psihat: complex d-dimensional array of Fourier coefficients
%   a:      left-hand points of the computational domain, in the form
%           a(1:end)
%   b:      right-hand points of the computational domain, in the form
%           a(1:end)
%   y:      cell array containing the column vectors y{:}, of the
%           coordinates of the points in the physical domain.
%           From y it possible to generate the grid points
%           [Y{1:end}] = ndgrid(y{1:end})
%
% Output variables:
%   psi:    complex nd-array of Fourier series evaluated at ndgrid(y{1:end})
N = size(psihat);
E = cell(1,2);
for d = 1:2
  E{d} = exp(1i*2*pi*(y{d}-a(d))/(b(d)-a(d))*(-N(d)/2:N(d)/2-1))/...
         sqrt(b(d)-a(d));
end
psi = E{1}*psihat*E{2}.';
%!test
%! a = [-20,-20];
%! b = [20,20];
%! N = [64,64];
%! M = N;
%! psihat = 2*rand(N)-1+1i*(2*rand(N)-1);
%! psihaty = ifft2(ifftshift(psihat))/prod(sqrt(b-a)./N);
%! ref = psihaty;
%! for d = 1:2
%!   y{d} = linspace(a(d),b(d),M(d)+1)';
%!   y{d} = y{d}(1:M(d));
%! end
%! psihaty = igridft2(psihat,a,b,y);
%! assert(psihaty(:),ref(:),1024*eps)
