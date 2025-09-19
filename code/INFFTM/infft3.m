function psi = infft3(psihat,a,b,Xi,loop)
% Function file: psi = infft3(psihat,a,b,Xi)
% Function file: psi = infft3(psihat,a,b,Xi,loop)
%
% Description: Evaluate 3d-array of Fourier coefficients by NFFT.
%
% In case NFFT is not available, a SLOW for-loop is used.
%
% Input variables:
%   psihat: complex 3d-array of Fourier coefficients
%   a:      row vector of left-hand points of the computational domain
%   b:      row vector of right-hand points of the computational domain
%   Xi:     2d-array containing in row d the d-th components of the
%           points, with Xi(d,:) in [a(d),b(d)), d=1:3
%   loop:   if true, the SLOW for-loop implementation (without NFFT)
%           is used (optional, default = false)
%
% Output variables:
%   psi:    complex row vector of the solution evaluated at Xi
nfftpath
if (nargin == 4)
  loop = false;
end
N = size(psihat);
M = size(Xi,2);
if (havenfft && ~loop)
  zeta = zeros(3,size(Xi,2));
  for d = 1:3
    zeta(d,:) = mod((Xi(d,:)-a(d))/(a(d)-b(d)),1)-0.5; % scale into (-1/2,1/2]
  end
  cs1 = ones(N(1),1);
  cs1(2-mod(N(1)/2,2):2:end) = -1;
  cs2 = ones(N(2),1);
  cs2(2-mod(N(2)/2,2):2:end) = -1;
  cs3 = ones(N(3),1);
  cs3(2-mod(N(3)/2,2):2:end) = -1;
  cs = kron(kron(cs1,cs2),cs3);
  plan = nfft_init_3d(N(1),N(2),N(3),M);
  % plan = nfft_init_guru(3,N(1),N(2),N(3),M,...
  %                       2^ceil(log2(2*N(1))),2^ceil(log2(2*N(2))),...
  %                       2^ceil(log2(2*N(3))),8,...
  %                       PRE_PSI,FFTW_ESTIMATE);
  psihat = permute(psihat,[3,2,1]);
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
  warning('Slow for-loop evaluation. NFFT installation recommended.')
  psi = complex(zeros(1,M));
  E = cell(1,3);
  for m = 1:M
    E{1}(:,1,1) = exp(1i * 2 * pi * (Xi(1,m) - a(1)) / (b(1) - a(1)) * ...
                      (-N(1)/2:N(1)/2 - 1)) / sqrt(b(1) - a(1));
    E{2}(1,:,1) = exp(1i * 2 * pi * (Xi(2,m) - a(2)) / (b(2) - a(2)) * ...
                      (-N(2)/2:N(2)/2 - 1)) / sqrt(b(2) - a(2));
    E{3}(1,1,:) = exp(1i * 2 * pi * (Xi(3,m) - a(3)) / (b(3) - a(3)) * ...
                      (-N(3)/2:N(3)/2 - 1)) / sqrt(b(3) - a(3));
    EE = bsxfun(@times,E{1} * E{2},E{3});
    psi(m) = sum(psihat(:) .* EE(:));
  end
end
%!demo
%! a = [-20,-20,-20];
%! b = [20,20,20];
%! N = [64,64,64];
%! psihat = 2*rand(N)-1+1i*(2*rand(N)-1);
%! tic
%! psihatXi = ifftn(ifftshift(psihat))/prod(sqrt(b-a)./N);
%! ref = psihatXi;
%! disp('ifft')
%! toc
%! M = prod(N);
%! for d = 1:3
%!   temp{d} = linspace(a(d),b(d),N(d)+1)';
%!   temp{d} = temp{d}(1:N(d));
%! end
%! [Temp{1:3}] = ndgrid(temp{1:3});
%! Xi(1,:) = Temp{1}(:)';
%! Xi(2,:) = Temp{2}(:)';
%! Xi(3,:) = Temp{3}(:)';
%! tic
%! psihatXi = infft3(psihat,a,b,Xi);
%! disp('NFFT')
%! error_inf = norm(psihatXi.'-ref(:),inf)
%! toc
%! disp('')
%! M = 1000;
%! Xi = rand(3,M);
%! for d = 1:3
%!   Xi(d,:) = a(d)+(b(d)-a(d))*Xi(d,:);
%! end
%! tic
%! psihatXi = infft3(psihat,a,b,Xi,true);
%! ref = psihatXi;
%! disp('one-loop')
%! toc
%! tic
%! psihatXi = infft3(psihat,a,b,Xi);
%! disp('NFFT')
%! error_inf = norm(psihatXi.'-ref(:),inf)
%! toc
%!test
%! a = [-20,-20,-20];
%! b = [20,20,20];
%! N = [16,16,16];
%! psihat = 2*rand(N)-1+1i*(2*rand(N)-1);
%! psihatXi = ifftn(ifftshift(psihat))/prod(sqrt(b-a)./N);
%! ref = psihatXi;
%! M = prod(N);
%! for d = 1:3
%!   temp{d} = linspace(a(d),b(d),N(d)+1)';
%!   temp{d} = temp{d}(1:N(d));
%! end
%! [Temp{1:3}] = ndgrid(temp{1:3});
%! Xi(1,:) = Temp{1}(:)';
%! Xi(2,:) = Temp{2}(:)';
%! Xi(3,:) = Temp{3}(:)';
%! psihatXi = infft3(psihat,a,b,Xi);
%! assert(psihatXi.',ref(:),256*eps)
