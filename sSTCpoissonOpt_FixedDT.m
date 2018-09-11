function [s S c] = sSTCpoissonOpt_FixedDT(D,T,Kr,K0,L,lamda,h,p,p2, S0)
% with fixed D=S-s, and T, the (s,S,T) cost is convex in S and so we
% compute the minimizing S. Then, s = S-D.
if nargin < 9
    p2 = 0;
end
if nargin < 10
    S0=1;
end
S = S0;  % initial estimate
s = S-D;
c = sSTCpoisson(s,S,T,Kr,K0,L,lamda,h,p,p2);
c2 = sSTCpoisson(s+1,S+1,T,Kr,K0,L,lamda,h,p,p2);
pos = false;
while c < c2
    pos = true;
    S = S-1;
    s = s-1;
    c2 = c;
    c = sSTCpoisson(s,S,T,Kr,K0,L,lamda,h,p,p2);
end
if pos
    s = s+1;
    S = S+1;
    c = c2;
else  % pos was false, so we need to increase s and S
    while c > c2
        S = S+1;
        s = s+1;
        c = c2;
        c2 = sSTCpoisson(s+1,S+1,T,Kr,K0,L,lamda,h,p,p2);
    end
    %s = s-1;
    %S = S+1;
    %c = c2;
end
end

