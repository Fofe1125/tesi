
function [x, h] = neqdeltagrid(a,b,delta,m)
% [x,h] = nqdeltagrid(a,b,delta,m)
% 
% Generate non equispaced vector in an interval [a,b] with
% x(1) = a, x(m) = b, h(i) = x(i+1)-x(i), h(i)=(1+delta)h(i-1)
%
x = NaN(1,m);
h = NaN(1,m-1);

x(1) = a;
if delta > 0
    h(1) = (delta*(b-a))/((1 + delta)^(m-1) - 1);
else
    h(1) = (b-a)/(m-1);
end

for i = 2:m-1

    h(i) = (1 + delta)*h(i-1);
    x(i) = x(i-1) + h(i-1);

end

x(m) = b;
h(m-1) = b - x(m-1);

end
