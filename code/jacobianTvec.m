function JTv = jacobianTvec(F, x, v)

n = numel(x);

JTv = zeros(n,1);
epsFD = 1e-8; 

Fx = F(x); 

for i = 1:n

    dx = zeros(n,1);
    dx(i) = epsFD;
    Fxp = F(x + dx);
    JTv(i) = v' * (Fxp - Fx) / epsFD;

end

end