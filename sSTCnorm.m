function [y l u p0]=sSTCnorm(s,S,T,Kr,K0,L,mi,sigma,h,p,p2,eps)
% p2 is the p in the Hadley-Whittin notation
% p is the phat in the Hadley-Whittin notation
%disp(['calling sSTCnorm(s=' num2str(s) ', S=' num2str(S) ', T=' num2str(T) ')']);  % itc: HERE rm asap
if nargin < 11
    p2 = 0;
end
% eps is the convergence accuracy
if nargin <12
    eps = 10^-8;
end
% convert into Hadley-Whittin notation
J=Kr;
A=K0;
r=s;
R=S;
l=mi;
m=mi*L;
t=L;
phat=p;
IC=h;

K1 = J/T;
nom = sum1(r,R,T,t,l,sigma,m,IC,phat,p2,eps) + H(R,T,t,l,sigma,m,IC,phat,p2);
denom = sum2(r,R,T,l,sigma,eps) + complnormcdf(R-r, l.*T, sigma.*sqrt(T));

y = K1 + (A+nom)/(T*denom);
l = K1 + nom/(T*denom);
u = K1 + A/T + nom/(T*denom);
p0 = 1/denom;
end


% auxiliary functions


function y = sum1(r,R,T,t,l,sigma,m,IC,phat,p2,eps)
y = 0;
n=1;
last = 0; count = 0;
while true
%    sum=0;
%    for j=1:R-r
%        sum = sum + poisspdf(R-r-j,n*l*T)*H(r+j,T,t,l,m,IC,phat,p2);
%    end
    sum = quadvec(@(x) finner2(x,r,R,T,l,sigma,t,n,m,IC,phat,p2), 0, R-r);
    y = y+sum;
    last = last+sum;
    count = count + 1;
    if count == 4
        if abs(last/y) < eps
            break;  % series converged
        else
            %disp(['sum1=' num2str(sum) ' y=' num2str(y) ' n=' num2str(n)]);
            count = 0;
            last = 0;
        end
    end
    n = n+1;
end
%disp(['sum1=' num2str(y) ' n=' num2str(n)]);
end


function y = H(rpj,T,t,l,sigma,m,IC,phat,p2)
z = 0;
if nargin == 9 && p2~=0
    z = p2.*eP(rpj,T,l,sigma,t);
end
y = IC.*T.*(rpj - m - l.*T./2) + (IC+phat).*bP(rpj,T,l,sigma,t) + z;
%disp(['H(' num2str(rpj) ')=' num2str(y)]);
end


function y = sum2(r,R,T,l,sigma,eps)
y = 0;
n=2;
last = 0; count = 0;
while true
%    sum=0;
%    for j=1:R-r
%        sum = sum + n*poisspdf(R-r-j,(n-1)*l*T)*complpoisscdf(j,l*T);
%    end
    sum = quadvec(@(x) finner3(x,r,R,T,l,sigma,n), 0, R-r);
    y = y+sum;
    last = last+sum;
    count = count + 1;
    if count == 4
        if abs(last/y) < eps
            break;  % series converged
        else
            %disp(['sum2=' num2str(sum) ' y=' num2str(y) ' n=' num2str(n)]);
            count = 0;
            last = 0;
        end
    end
    n = n+1;
end
%disp(['sum2=' num2str(y) ' n=' num2str(n)]);
end


function y=complnormcdf(x,mt,st)
y = 1.0 - normcdf(x,mt,st);
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
%y = quadgk(@(x) finta(x, rpj, l, sigma), t, t+T);
x = rpj; mi=l; L=t;
t1 = (sigma.^2).*(W1(x,L+T,mi,sigma) - W1(x,L,mi,sigma));
t2 = mi.*(V1(x,L+T,mi,sigma)-V1(x,L,mi,sigma)) - x.*(V0(x,L+T,mi,sigma)-V0(x,L,mi,sigma));
y = t1 + t2;
end

function y=V0(x,T,m,s)
xstn = (x-m.*T)./(s.*sqrt(T));
xstp = (x+m.*T)./(s.*sqrt(T));
t1 = (T-x./m-(s.^2)/(2.*m.^2)).*complnormcdf(xstn,0,1);
t2 = s.*sqrt(T)*normpdf(xstn,0,1)./m;
t3 = (s.^2).*exp(2.*m.*x./(s.^2)).*complnormcdf(xstp,0,1)./(2.*m.^2);
y = t1+t2+t3;
end

function y=V1(x,T,m,s)
xstn = (x-m.*T)./(s.*sqrt(T));
xstp = (x+m.*T)./(s.*sqrt(T));
t1 = complnormcdf(xstn,0,1).*(T.^2 - (x./m).^2 - 2.*(s.^2).*x./(m.^3) - 3.*(s.^4)./(2.*m.^4))./2;
t2 = normpdf(xstn,0,1).*s.*sqrt(T).*(m.*T + 3.*s.*s./m + x)./(2.*m.^2);
t3 = (s.^2).*(x - 3.*(s.^2)./(2.*m)).*exp(2.*m.*x./(s.^2)).*complnormcdf(xstp,0,1)./(2.*m.^3);
y = t1 + t2 - t3;
end

function y=W1(x,T,m,s)
xstn = (x-m.*T)./(s.*sqrt(T));
xstp = (x+m.*T)./(s.*sqrt(T));
t1 = (s.^2).*(1 + m.*x./(s.^2)).*complnormcdf(xstn,0,1)./(m.^3);
t2 = 2.*s.*sqrt(T).*normpdf(xstn,0,1)./(m.^2);
t3 = (x-(s.^2)./m).*exp(2.*m.*x./(s.^2)).*complnormcdf(xstp,0,1)./(m.^2);
y = t1 - t2 + t3;
end


function y = finta(ksi,rpj,mi,sigma)
quant = (rpj-mi.*ksi)./(sigma.*sqrt(ksi));
y1 = sigma.^2.*ksi.*normpdf(quant);
y2 = (mi.*ksi-rpj).*sigma.*sqrt(ksi).*complnormcdf(quant,0,1);
y = (y1+y2)./(sigma.*sqrt(ksi));
end


function y = finner(t,x,rpj,mi,sigma)
%disp(['t=' num2str(t) ' x=' num2str(x)]);
nom1 = t-mi.*x;
denom1 = sigma.*sqrt(x);
%disp(['t-mi.*x=' num2str(nom1)]);
%disp(['sigma.*sqrt(x)=' num2str(denom1)]);
varg = nom1./denom1;
v = normpdf(varg);
denom = sigma.*sqrt(x);
ft = t-rpj;
%disp(['varg=' num2str(varg) ' v=' num2str(v) ' ft=' num2str(ft) ' denom=' num2str(denom)]);
y = (ft.*v)./denom;
%disp(['finner(' num2str(t) ',' num2str(x) ',...)=' num2str(y)]);
end


function y = fint(x,rpj,mi,sigma,ulimit)
y = 0;
z = 1;
a = rpj;
b = ulimit;
while abs(z) > 10^-12
    z = quadvec(@(t) finner(t,x,rpj,mi,sigma), a, b);
    %disp(['a=' num2str(a) ' b=' num2str(b) ' y=' num2str(y) ' z=' num2str(z)]);
    y = y+z;
    a = b;
    b = b + ulimit;
end
end


function y=finner2(x,r,R,T,l,sigma,t,n,m,IC,phat,p2)
y = H(r+x,T,t,l,sigma,m,IC,phat,p2).*normpdf(R-r-x, n.*l.*T, sigma.*sqrt(n.*T));
end

function y = finner3(x,r,R,T,l,sigma,n)
y = n.*complnormcdf(x, l.*T, sigma.*sqrt(T)).*normpdf(R-r-x, (n-1).*l.*T, sigma.*sqrt((n-1).*T));
end

