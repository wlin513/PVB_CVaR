function  y = exponential_utility(x, k)

m = 0.5;
s = 0.5;

if k == 0
    f = @(x,k) x;
else
    f = @(x,k) (1 - exp(-k * x)) / k;
end

y = (x - m)/s;
ndx = not(y == 0);
y(ndx) = f(y(ndx), k);

y = y*s + m;

end