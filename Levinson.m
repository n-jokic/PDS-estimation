function [a, s] = Levinson(x, p)
N = length(x);
Rx = xcorr(x)/length(x);
    off = length(Rx) - 1;
    off = off/2 + 1;
a = zeros(1, p);
a(1) = -Rx(1+off)/Rx(off);
P = (1-abs((a(1)))^2)*Rx(off);
anew = zeros(1, p);
for k = 2 : p


    a(k) = -(Rx(k+ off) + sum(a(1: k -1).*Rx(k + off - 1:-1:1+ off)))/P;
    
    
    for i = 1 : k-1
        anew(i) = a(i) + a(k)*conj(a(k - i));
    end
    anew(k) = a(k);
    a = anew;
    P = (1- abs(a(k))^2)*P;
end


s  = Rx(off) + a*Rx(off - 1 :-1: off-p)';
a = [1 a];
end


