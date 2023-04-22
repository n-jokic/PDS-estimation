function [a, s] = covar(x, p)
N = length(x);
cxx = zeros(p+1, p+1);

for j = 1 : p + 1
    for k = 1 : p + 1
        cxx(j, k) = conj(x(p - j + 1 + 1 : N - j + 1))*(x(p - k + 1 + 1 : N - k + 1))'/(N-p);
    end
end


a = -cxx(2:p+1, 2:p+1)\cxx(2:p+1, 1);
a = a';
s = cxx(1, 1) + sum(a.*cxx(1, 2:p+1));
a = [1, a];
end

