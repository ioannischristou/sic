function [s Q T c] = snQTCnormOptNonLin(Kr,K0,L,mi,sigma,h,p)
T0=(3*sigma/mi)^2;
r0=(mi+3*sigma)*(L+T0);
Q0=10;
v0 = [r0 Q0 T0];
[v c]= fmincon(@(x) snQTCnormV(x,Kr,K0,L,mi,sigma,h,p),v0,[],[],[],[],[-Inf 1e-2 T0],[Inf Inf Inf]);
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