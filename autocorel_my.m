function [a, s] = autocorel_my(x, p)
N = length(x);
Rx = zeros(2*p+1, 1);

for k = (p+1) : (2*p+1)
    Rx(k) = conj(x(1:N -(k - p - 1)))*(x((k - p - 1) + 1 : N))'/N;
end

Rx(1:p) = conj(Rx(end : -1 : p+2));

R = toeplitz(Rx(p+1 : end  - 1), Rx(p+1: -1 : 2));

a = -R^-1 * Rx(p+2 : end);
s = Rx(p+1) + Rx(p : -1 : 1)'*a;
a = [1 a'];

end
