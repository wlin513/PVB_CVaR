function  y = power_utility(x, k)

m = 0.5;
s = 0.5;

f = @(x, k) x./abs(x) .* abs(x).^k; % exp compression

y = (x - m)/s;
ndx = not(y == 0);
y(ndx) = f(y(ndx), k);

y = y*s + m;

end