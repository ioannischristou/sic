function y  = sSTCnorm_2(v,Kr,K0,L,mi,sigma,h,p,p2,eps)
if nargin < 9
    p2=0;
end
if nargin < 10
    eps=10^-8;
end
epsq = 10^-3;
tnot = (3*sigma/mi)^2;
s = v(1,1);
S = v(2,1);
T = v(3,1);
% bring into the feasible region
if T<tnot
    T=tnot;
end
if S-s<=0
    s = S-epsq;
end
y = sSTCnorm(s,S,T,Kr,K0,L,mi,sigma,h,p,p2,eps);
end

