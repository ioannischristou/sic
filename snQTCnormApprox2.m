function y = snQTCnormApprox2(s,Q,T,Kr,K0,L,mi,sigma,h,p,p2,Qmin)
if nargin < 12
    Qmin = 1.e-6;
end
if nargin < 11
    p2 = 0;
end
y = snQTCnorm(s,Q,T,Kr,0,L,mi,sigma,h,p,p2,Qmin);
P0 = mi*T/Q;
if P0 > 1
    P0 = 1;
end
y = y+K0.*P0./T;
end