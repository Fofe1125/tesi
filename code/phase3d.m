function theta = phase3d(X,Y,Z,centerLine,kappa,x0)
% Calcola la fase di un vortice 3D tramite Biot-Savart
% X,Y,Z: griglia 3D (stesse dimensioni)
% centerLine: Nx3 array con punti della linea vorticosa
% kappa: circolazione (default 2*pi)
% x0: punto di riferimento (default [0,0,0])

    if nargin < 5
        kappa = 2*pi;
    end
    if nargin < 6
        x0 = [0,0,0];
    end

    ds = diff(centerLine,1,1);
    s_mid = (centerLine(1:end-1,:) + centerLine(2:end,:))/2; 

    [nx,ny,nz] = size(X);
    theta = zeros(nx,ny,nz);

    lambda = linspace(0,1,5);
    w = [1 4 6 4 1]/16;

    parfor i = 1:nx
        for j = 1:ny
            for k = 1:nz
                x = [X(i,j,k), Y(i,j,k), Z(i,j,k)];
                dx = x - x0;
                Sval = 0;
                for q = 1:length(lambda)
                    xp = x0 + lambda(q)*dx;
                    u = vField(xp,s_mid,ds,kappa);
                    Sval = Sval + w(q)*dot(u,dx);
                end
                theta(i,j,k) = Sval;
            end
        end
    end
end

function u = vField(x,s_mid,ds,kappa)
% Biot-Savart velocity at point x
    u = [0,0,0];
    parfor j = 1:size(ds,1)
        r = x - s_mid(j,:);
        rnorm = norm(r);
        if rnorm > 1e-6
            u = u + cross(ds(j,:),r)/(rnorm^3);
        end
    end
    u = (kappa/(4*pi)) * u;
end
