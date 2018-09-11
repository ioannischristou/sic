function [s Q T c] = snQTCnormOptAD(Kr,K0,L,mi,sigma,h,p)
T0=(3*sigma/mi)^2;
r0=(mi+3*sigma)*(L+T0);
Q0=10;
v0 = [r0 Q0 T0];
epsv = [0.01 0.1 0.01];
%tryorder = [1 2 1 3];
tryorder = [1 2 3];
[v c]= nonlinoptAD(@(x) snQTCnormV(x,Kr,K0,L,mi,sigma,h,p),v0, tryorder, epsv, 100);
s = v(1);
Q = v(2);
T = v(3);
c = snQTCnorm2(s,Q,T,Kr,K0,L,mi,sigma,h,p);
end

function y=snQTCnormV(v,Kr,K0,L,mi,sigma,h,p)
s = v(1);
Q=v(2);
T=v(3);
y = snQTCnorm2(s,Q,T,Kr,K0,L,mi,sigma,h,p);
end