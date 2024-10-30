function [a,b] = crack_integer(m,n)
integer = m*n;
start = floor(sqrt(integer));
factor = integer/start;
while factor ~= floor(factor)
    start = start - 1;
    factor = integer/start;
end
if m <= n
    a = start;
    b = factor;
else
    a = factor;
    b = start;
end
end

