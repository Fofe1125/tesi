function [x, h, m] = neqdeltagrid0(xc, xL, delta, h1)

% compute m
m = ceil(log(delta*xc/h1 + 1)/log(1 + delta) + 1);

% preallocate for code speed
x = NaN(1,3*m);
h = NaN(1,3*m);

% central point 
x(m) = xc;
h(m) = h1;

% grid to 0
for i = m-1:-1:2
    h(i) = (1 + delta)*h(i+1);
    x(i) = x(i+1) - h(i);
end

% from xc to xL
i = m + 1;
while x(i - 1) < xL
    h(i) = (1 + delta)*h(i-1);
    x(i) = x(i-1) + h(i);
    i = i + 1;
end

x = x(~isnan(x));
h = h(~isnan(h));

% only positive values
x = x(x >= 0);

% always add 0
x = [0, x];
h = diff(x);

% check for the nearest point to 0
if numel(x) > 2
    if h(1) < 0.7*h(2)
        x(2) = [];
        h = diff(x);
    end
end

m = numel(x);
end
