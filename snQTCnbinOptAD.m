function [s Q T c] = snQTCnbinOptAD(Kr,K0,L,lamda,pl,h,p)
T0=0.1;
r0=1;
Q0=10;
v0 = [r0 Q0 T0];
epsv = [1 1 0.01];
%tryorder = [1 2 1 3];
tryorder = [1 2 3];
[v c]= nonlinoptAD(@(x) snQTCnbinV(x,Kr,K0,L,lamda,pl,h,p),v0, tryorder, epsv, 100);
s = v(1);
Q = v(2);
T = v(3);
c = snQTCnbin(s,Q,T,Kr,K0,L,lamda,pl,h,p);
end

function y=snQTCnbinV(v,Kr,K0,L,lamda,pl,h,p)
s = v(1);
Q=v(2);
if Q<=0
    Q=1;
end
T=v(3);
if T<=0
    T=0.01;
end
y = snQTCnbin(s,Q,T,Kr,K0,L,lamda,pl,h,p);
end