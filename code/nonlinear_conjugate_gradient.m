function [x, its] = nonlinear_conjugate_gradient(F,JF,x0,tol)

    d0 = -JF(x0);
    g = d0;
    gold = g;
    dold = d0;
    xold = x0;
    x = xold;
    its = 0;

    while norm(dold, 2) >= tol*norm(d0, 2) 

        its = its + 1;
        
        f = @(alpha) F(xold + alpha*dold);
        alfa = golden_search(f,0,1e-2,1e-12);
        x = xold + alfa*dold;
        g = -JF(x);
        beta = (g.'*JF(x))/(gold.'*JF(xold));
        d = g + beta*dold;

        dold = d;
        xold = x;
        gold = g;

    end
end
