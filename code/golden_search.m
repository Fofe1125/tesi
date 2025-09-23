function xmin = goldenSearch(f, a, b, tol)
    % goldenSearch - Minimizzazione con Golden Section Search
    %
    % Sintassi:
    %   xmin = goldenSearch(f, a, b, tol)
    %
    % Input:
    %   f   - handle della funzione (es: @(x) (x-2).^2 + 1)
    %   a   - estremo sinistro dell'intervallo
    %   b   - estremo destro dell'intervallo
    %   tol - tolleranza (default 1e-6)
    %
    % Output:
    %   xmin - stima del punto di minimo

    if nargin < 4
        tol = 1e-6;
    end

    phi = (sqrt(5) - 1) / 2;   % rapporto aureo ≈ 0.618
    x1 = b - phi * (b - a);
    x2 = a + phi * (b - a);
    f1 = f(x1);
    f2 = f(x2);

    while (b - a) > tol
        if f1 > f2
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + phi * (b - a);
            f2 = f(x2);
        else
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = b - phi * (b - a);
            f1 = f(x1);
        end
    end

    xmin = (a + b) / 2;
end