function v = computeKronProd(A,B,u)

    m = size(A,1);
    n = size(B,2);

    P = reshape(u, m,n);

    V = B*P*A.';

    v = reshape(V,m*n,1);

end
