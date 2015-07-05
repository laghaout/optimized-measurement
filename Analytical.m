function [X, Y, P] = Analytical(alpha, t, x, y, n)

if alpha == 0
    alpha = 1e-12;
    warning('Setting alpha = %f instead of alpha = 0 to avoid underflow\n', alpha);
end

if n > 0
    P = 0;
    for m = 1:10*n
        P = P + (exp(-abs(alpha)^2)*abs(alpha)^(2*m)/factorial(m))*(abs(x)^2 + abs(y)^2*(1-t)*abs(m/alpha - conj(alpha))^2 + conj(x)*y*sqrt(1-t)*(m/alpha - conj(alpha)) + x*conj(y)*sqrt(1-t)*(m/conj(alpha) - alpha) + abs(y)^2*t);
    end
else
    m = n;
    P = (exp(-abs(alpha)^2)*abs(alpha)^(2*m)/factorial(m))*(abs(x)^2 + abs(y)^2*(1-t)*abs(m/alpha - conj(alpha))^2 + conj(x)*y*sqrt(1-t)*(m/alpha - conj(alpha)) + x*conj(y)*sqrt(1-t)*(m/conj(alpha) - alpha) + abs(y)^2*t);
end

OverallFactor = (1/sqrt(P))*exp(-abs(alpha)^2/2)*alpha^n/sqrt(factorial(n));

X = OverallFactor*(x + y*sqrt(1-t)*(n/alpha - conj(alpha)));

Y = OverallFactor*(y*sqrt(t));
