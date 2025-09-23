function JF = approximate_jacobian(F, x0, eps)
%
%   Compute the action of the Jacobian matrix F evaluated 
%   at x0 on the vector u
%

if nargin < 3
    eps = 1e-6;
end

n = length(x0);
JF = speye(length(x0));
I = speye(n);

for i = 1:n

    y = I(:,i);
    JF(:,i) = imag(F(x0 + 1i*eps*y)/eps);

end

if n < 50
    JF = full(JF);
end

end