function [s S c] = sSTCnormOpt_FixedT(T,Kr,K0,L,mi,sigma,h,p, p2, epsd)
if nargin < 10
    epsd = 1;
end
if nargin < 9
    p2 = 0;
end
% step 1
ystar = findGmin(T,L,mi,sigma,h,p,p2,epsd);
s = ystar;
S_0 = ystar;
while true
    s = s-epsd;
    c = sSTCnorm(s,S_0,T,0,K0,L,mi,sigma,h,p,p2);  % set Kr=0 to work
    Gs = G(s,T,L,mi,sigma,h,p,p2);
    if c<=Gs
        disp(['s=' num2str(s) ' c=' num2str(c) ' G(s)=' num2str(Gs)]);
        break;
    end
end
s0 = s;
c0 = c;
Su0 = S_0;
S = Su0 + epsd;
GS = G(S,T,L,mi,sigma,h,p,p2);
disp(['c0=' num2str(c0) ' G(S)=' num2str(GS)]);
% step 2
while GS <= c0
    css = sSTCnorm(s,S,T,0,K0,L,mi,sigma,h,p,p2);  % set Kr=0 to work
    if css < c0
        Su0 = S;
        csS0 = sSTCnorm(s,Su0,T,0,K0,L,mi,sigma,h,p,p2);  % Kr=0 to work
        Gsp1 = G(s+epsd,T,L,mi,sigma,h,p,p2);
        while csS0 <= Gsp1
            s = s+epsd;
            csS0 = sSTCnorm(s,Su0,T,0,K0,L,mi,sigma,h,p,p2);  % Kr=0 to work
            Gsp1 = G(s+epsd,T,L,mi,sigma,h,p,p2);  % this loop can be enhanced
        end
        c0 = csS0;
    end
    S = S+epsd;
    GS = G(S,T,L,mi,sigma,h,p,p2);
    disp(['s=' num2str(s) ' S=' num2str(S) ' c0=' num2str(c0) ' G(S)=' num2str(GS)]);
end
S = Su0;
c = sSTCnorm(s,S,T,Kr,K0,L,mi,sigma,h,p,p2);  % report real value
end

function z = findGmin(T,L,mi,sigma,h,p,p2,epsd)
y = 0;
c = G(y,T,L,mi,sigma,h,p,p2);
c2 = G(y+epsd,T,L,mi,sigma,h,p,p2);
if c==c2
    z=y;
    return;
end
if c<c2
    z = y;
    while true
        y = y-epsd;
        c2 = G(y,T,L,mi,sigma,h,p,p2);
        if c2>=c
            z = y+epsd;
            break;
        else
            c = c2;
        end
    end
else
    z = y;
    while true
        y = y+epsd;
        c2 = G(y,T,L,mi,sigma,h,p,p2);
        if c2>=c
            z = y-epsd;
            break;
        else
            c = c2;
        end
    end
end
end

function y = G(s,T,L,mi,sigma,h,p,p2)
z = 0;
if nargin == 8 && p2~=0
    z = p2*eP(s,T,mi,sigma,L);
end
y = h*T*(s - L*mi - mi*T/2) + (h+p)*bP(s,T,mi,sigma,L) + z;
y = y/T;  % need to divide by T to account for the review period length
end

function y = eP(rpj,T,l,sigma,t)
fh = @(x) (x-rpj).*(normpdf(x,l.*(t+T),sigma.*sqrt(t+T))-normpdf(x,l.*T,sigma.*sqrt(T)));
%y = quadl(fh, rpj, +Inf);
stepsize = 100;
a = rpj;
b = rpj+stepsize;
z=1;
y=0;
while abs(z)>10^-12
   z = quadvec(fh, a, b);
   y = y+z;
   a = b;
   b = b+stepsize;
end
disp(['eP(' num2str(rpj) ')=' num2str(y)])
end

function y = bP(rpj,T,l,sigma,t)
y = quadvec(@(x) finta(x, rpj, l, sigma), t, t+T);
end

function y = finta(ksi,rpj,mi,sigma)
quant = (rpj-mi.*ksi)./(sigma.*sqrt(ksi));
y1 = sigma.^2.*ksi.*normpdf(quant);
y2 = (mi.*ksi-rpj).*sigma.*sqrt(ksi).*complnormcdf(quant,0,1);
y = (y1+y2)./(sigma.*sqrt(ksi));
end

%function y = bP(rpj,T,l,t)
%t1  = l*(((t+T)^2)*complpoisscdf(rpj-1,l*(t+T))-(t^2)*complpoisscdf(rpj-1,l*t))/2;
%t2 = (complpoisscdf(rpj+1,l*(t+T))-complpoisscdf(rpj+1,l*t))*rpj*(rpj+1)/(2*l);
%t3 = rpj*((t+T)*complpoisscdf(rpj,l*(t+T))-t*complpoisscdf(rpj,l*t));
%y = t1 + t2 - t3;
%end

%function y = eP(rpj,T,l,t)
%t1 = l*(t+T)*complpoisscdf(rpj-1,l*(t+T)) - rpj*complpoisscdf(rpj,l*(t+T));
%t2 = l*t*complpoisscdf(rpj-1,l*t) + rpj*complpoisscdf(rpj,l*t);
%y = t1 - t2;
%end

function y=complnormcdf(x,mt,st)
y = 1.0 - normcdf(x,mt,st);
end