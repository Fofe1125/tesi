function psi = igridftn(psihat,a,b,y)
% Function file: psi = igridftn(psihat,a,b,y)
%
% Description: Evaluate nd-array of Fourier coeffs to a rectilinear grid.
%
% Input variables:
%   psihat: complex nd-array of Fourier coefficients
%   a:      row vector of left-hand points of the computational domain
%   b:      row vector of right-hand points of the computational domain
%   y:      cell array containing the column vectors y{:}, of the
%           coordinates of the points in the physical domain.
%           From y it possible to generate the grid points
%           [Y{1:end}] = ndgrid(y{1:end})
%
% Output variables:
%   psi:    complex nd-array of Fourier series evaluated at ndgrid(y{1:end})
N = size(psihat);
dd = sum(N>1);
E = cell(1,dd);
for d = 1:dd
  E{d} = exp(1i*2*pi*(y{d}-a(d))/(b(d)-a(d))*(-N(d)/2:N(d)/2-1))/...
         sqrt(b(d)-a(d));
end
psi = ndcovlt(psihat,E);
%!demo
%! a = [-20,-20,-20];
%! b = [20,20,20];
%! N = [64,64,64];
%! M = N;
%! psihat = 2*rand(N)-1+1i*(2*rand(N)-1);
%! tic
%! psihaty = ifftn(ifftshift(psihat))/prod(sqrt(b-a)./N);
%! ref = psihaty;
%! disp('ifft')
%! toc
%! for d = 1:3
%!   y{d} = linspace(a(d),b(d),M(d)+1)';
%!   y{d} = y{d}(1:M(d));
%! end
%! tic
%! psihaty = igridftn(psihat,a,b,y);
%! disp('ndcovlt')
%! error_inf = norm(psihaty(:)-ref(:),inf)
%! toc
%! for d = 1:3
%!   E{d} = exp(1i*2*pi*(y{d}-a(d))/(b(d)-a(d))*(-N(d)/2:N(d)/2-1))/...
%!          sqrt(b(d)-a(d));
%! end
%! tic
%! psihaty = zeros(M);
%! for k3 = 1:N(3)
%!   temp = E{1}*(psihat(:,:,k3)*E{2}.');
%!   for m3 = 1:M(3)
%!     psihaty(:,:,m3) = psihaty(:,:,m3)+...
%!                       temp*E{3}(m3,k3);
%!   end
%! end
%! disp('two loops')
%! error_inf = norm(psihaty(:)-ref(:),inf)
%! toc
%!test
%! a = [-20,-20,-20];
%! b = [20,20,20];
%! N = [64,64,64];
%! M = N;
%! psihat = 2*rand(N)-1+1i*(2*rand(N)-1);
%! psihaty = ifftn(ifftshift(psihat))/prod(sqrt(b-a)./N);
%! ref = psihaty;
%! for d = 1:3
%!   y{d} = linspace(a(d),b(d),M(d)+1)';
%!   y{d} = y{d}(1:M(d));
%! end
%! psihaty = igridftn(psihat,a,b,y);
%! assert(psihaty(:),ref(:),1024*eps)
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
%! psihaty = igridftn(psihat,a,b,y);
%! assert(psihaty(:),ref(:),1024*eps)
