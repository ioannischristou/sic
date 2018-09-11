function [s S c] = sSTCpoissonOpt_FixedT(T,Kr,K0,L,lamda,h,p,p2)
if nargin < 8
    p2 = 0;
end
% step 1
ystar = findGmin(T,L,lamda,h,p,p2);
disp(['for T=' num2str(T) ' ystar=' num2str(ystar)]);  % itc: HERE rm asap
s = ystar;
S_0 = ystar;
while true
    s = s-1;
    c = sSTCpoisson(s,S_0,T,0,K0,L,lamda,h,p,p2);  % set Kr=0 to work
    Gs = G(s,T,L,lamda,h,p,p2);
        disp(['s=' num2str(s) ' c=' num2str(c) ' G(s)=' num2str(Gs)]);  % itc: HERE rm asap
    if c<=Gs
        disp(['s=' num2str(s) ' c=' num2str(c) ' G(s)=' num2str(Gs)]);
        break;
    end
end
s0 = s;
c0 = c;
Su0 = S_0;
S = Su0 + 1;
GS = G(S,T,L,lamda,h,p,p2);
disp(['c0=' num2str(c0) ' G(S)=' num2str(GS)]);
% step 2
while GS <= c0
    css = sSTCpoisson(s,S,T,0,K0,L,lamda,h,p,p2);  % set Kr=0 to work
    if css < c0
        Su0 = S;
        csS0 = sSTCpoisson(s,Su0,T,0,K0,L,lamda,h,p,p2);  % Kr=0 to work
        Gsp1 = G(s+1,T,L,lamda,h,p,p2);
        while csS0 <= Gsp1
            s = s+1;
            csS0 = sSTCpoisson(s,Su0,T,0,K0,L,lamda,h,p,p2);  % Kr=0 to work
            Gsp1 = G(s+1,T,L,lamda,h,p,p2);  % this loop can be enhanced
        end
        c0 = csS0;
    end
    S = S+1;
    GS = G(S,T,L,lamda,h,p,p2);
    disp(['s=' num2str(s) ' S=' num2str(S) ' c0=' num2str(c0) ' G(S)=' num2str(GS)]);
end
S = Su0;
c = sSTCpoisson(s,S,T,Kr,K0,L,lamda,h,p,p2);  % report real value
end

function z = findGmin(T,L,lamda,h,p,p2)
tol = 10^-6;
y = 0;
c = G(y,T,L,lamda,h,p,p2);
c2 = G(y+1,T,L,lamda,h,p,p2);
%if c==c2  
%    z=y;
%    return;
%end
if c<c2-tol  % used to be c < c2
    z = y;
    while true
        y = y-1;
        c2 = G(y,T,L,lamda,h,p,p2);
        if c2>=c+tol  % used to be c2>=c
            z = y+1;
            break;
        else
            c = c2;
        end
    end
else
    z = y;
    while true
        y = y+1;
        c2 = G(y,T,L,lamda,h,p,p2);
        if c2>=c+tol  % used to be c2>=c
            z = y-1;
            break;
        else
            c = c2;
        end
    end
end
end

function y = G(s,T,L,lamda,h,p,p2)
z = 0;
if nargin == 7 && p2~=0
    z = p2*eP(s,T,lamda,L);
end
y = h*T*(s - L*lamda - lamda*T/2) + (h+p)*bP(s,T,lamda,L) + z;
y = y/T;  % need to divide by T to account for the review period length
end

function y = bP(rpj,T,l,t)
t1  = l*(((t+T)^2)*complpoisscdf(rpj-1,l*(t+T))-(t^2)*complpoisscdf(rpj-1,l*t))/2;
t2 = (complpoisscdf(rpj+1,l*(t+T))-complpoisscdf(rpj+1,l*t))*rpj*(rpj+1)/(2*l);
t3 = rpj*((t+T)*complpoisscdf(rpj,l*(t+T))-t*complpoisscdf(rpj,l*t));
y = t1 + t2 - t3;
end

function y = eP(rpj,T,l,t)
t1 = l*(t+T)*complpoisscdf(rpj-1,l*(t+T)) - rpj*complpoisscdf(rpj,l*(t+T));
t2 = l*t*complpoisscdf(rpj-1,l*t) - rpj*complpoisscdf(rpj,l*t);
y = t1 - t2;
end

function y=complpoisscdf(x,lm)
y = 1.0 - poisscdf(x-1,lm);
end