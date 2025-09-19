function psi = infft(psihat,a,b,Xi,matvec)
% Function file: psi = infft(psihat,a,b,Xi)
% Function file: psi = infft(psihat,a,b,Xi,matvec)
%
% Description: Evaluate 1d-array of Fourier coefficients by NFFT.
%
% In case NFFT is not available, a matrix-vector approach is used.
%
% Input variables:
%   psihat: complex column vector of Fourier coefficients
%   a:      left-hand point of the computational domain
%   b:      right-hand point
%   Xi:     row vector containing the points, with Xi(1,:) in [a,b)
%   matvec: if true, the matrix-vector implementation (without NFFT)
%           is used (optional, default = false)
%
% Output variables:
%   psi:    complex row vector of the solution evaluated at Xi
nfftpath
if (nargin == 4)
  matvec = false;
end
N = length(psihat);
M = size(Xi,2);
if (havenfft && ~matvec)
  zeta = mod((Xi-a)/(a-b),1)-0.5; % scaling into [-1/2,1/2)
  cs = ones(N,1);
  cs(2-mod(N/2,2):2:end) = -1;
  plan = nfft_init_1d(N,M);
  psihat = psihat(:);
  nfft_set_f_hat(plan,psihat / sqrt(prod(b-a)) .* cs);
  % node-dependent precomputation
  nfft_set_x(plan,zeta);
  nfft_precompute_psi(plan);
  % transform
  nfft_trafo(plan);
  % function values
  psi = nfft_get_f(plan).';
  nfft_finalize(plan);
else
  warning('Matrix-vector product.')
  E = exp(1i * 2 * pi * (Xi(:) - a) / (b - a) * (-N/2:N/2 - 1)) / ...
       sqrt(b - a);
  psi = E * psihat(:);
  psi = psi(:).';
end
%!demo
%! a = -20;
%! b = 20;
%! N = 1024;
%! psihat = 2*rand(N,1)-1+1i*(2*rand(N,1)-1);
%! tic
%! psihatXi = ifft(ifftshift(psihat))/prod(sqrt(b-a)./N);
%! ref = psihatXi;
%! disp('ifft')
%! toc
%! M = N;
%! Xi = linspace(a,b,N+1);
%! Xi = Xi(1:N);
%! tic
%! psihatXi = infft(psihat,a,b,Xi);
%! disp('NFFT')
%! error_inf = norm(psihatXi.'-ref,inf)
%! toc
%! disp('')
%! M = 1000;
%! Xi = rand(1,M);
%! Xi = a+(b-a)*Xi;
%! tic
%! psihatXi = infft(psihat,a,b,Xi,true);
%! ref = psihatXi;
%! disp('matrix-vector')
%! toc
%! tic
%! psihatXi = infft(psihat,a,b,Xi);
%! disp('NFFT')
%! error_inf = norm(psihatXi-ref,inf)
%! toc
%!test
%! a = -20;
%! b = 20;
%! N = 256;
%! psihat = 2*rand(N,1)-1+1i*(2*rand(N,1)-1);
%! psihatXi = ifft(ifftshift(psihat))/prod(sqrt(b-a)./N);
%! ref = psihatXi;
%! M = N;
%! Xi(1,1:N+1) = linspace(a,b,N+1)';
%! Xi = Xi(1:N);
%! psihatXi = infft(psihat,a,b,Xi);
%! assert(psihatXi.',ref,256*eps)
