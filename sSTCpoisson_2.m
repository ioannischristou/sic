function y = sSTCpoisson_2(v, Kr, Ko, L, lamda, h, p, p2)
s = round(v(1,1));
S = round(v(2,1));
T = v(3,1);
% bring into the feasible region
if T<5/lamda
    T=5/lamda;
end
if S-s<=0
    s=S-1;
end
y = sSTCpoisson(s, S, T, Kr, Ko, L, lamda, h, p, p2);
end

