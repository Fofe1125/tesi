function [x, h, m] = neqdeltagrid0(xc, xL, delta, h1)
% neqdeltagrid0: costruisce una griglia non equispaziata
% che include sempre lo 0 e controlla il passo vicino a 0.

m = ceil(log(delta*xc/h1 + 1)/log(1 + delta) + 1);

x = NaN(1,3*m);
h = NaN(1,3*m);

% punto centrale
x(m) = xc;
h(m) = h1;

% griglia verso sinistra
for i = m-1:-1:2
    h(i) = (1 + delta)*h(i+1);
    x(i) = x(i+1) - h(i);
end

% griglia verso destra
i = m + 1;
while x(i - 1) < xL
    h(i) = (1 + delta)*h(i-1);
    x(i) = x(i-1) + h(i);
    i = i + 1;
end

% ripulisci NaN
x = x(~isnan(x));
h = h(~isnan(h));

% prendi solo parte positiva
x = x(x >= 0);

% aggiungi sempre lo 0 come primo punto
x = [0, x];
h = diff(x);

% --- controllo passo vicino a 0 ---
if numel(x) > 2
    if h(1) < 0.7*h(2)
        % elimina il primo punto positivo
        x(2) = [];
        h = diff(x);
    end
end

m = numel(x);
end
