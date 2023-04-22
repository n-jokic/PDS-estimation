function [a, s] = autocorel(x, p)
N = length(x);
Rx = zeros(2*p+1, 1);

for k = (p+1) : (2*p+1)
    Rx(k) = sum(conj(x(1:N-k+p+1)) .* x(1+k-(p+1):N))/N; 
end

Rx(p:-1:1) = conj(Rx(p+2:end));

R = toeplitz(Rx(p+1:2*p), Rx((p+1):-1:2));

a = -R^-1*Rx(p+2:end);
s  = Rx(p+1) + sum(a.*Rx(p:-1:1));
a = [1 transpose(a)];
end

