function [y l u p0]=sSTCnorm(s,S,T,Kr,K0,L,mi,sigma,h,p,p2,eps)
% p2 is the p in the Hadley-Whittin notation
% p is the phat in the Hadley-Whittin notation
disp(['calling sSTCnorm(s=' num2str(s) ', S=' num2str(S) ', T=' num2str(T) ';L=' num2str(L) ' mi=' num2str(mi) ' sigma=' num2str(sigma) ' h=' num2str(h) ' p=' num2str(p) ')']);  % itc: HERE rm asap
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
sum2v = sum2(r,R,T,l,sigma,eps);
cmpl = complnormcdf(R-r, l.*T, sigma.*sqrt(T));
denom = sum2v + cmpl;

y = K1 + (A+nom)/(T*denom);
l = K1 + nom/(T*denom);
u = K1 + A/T + nom/(T*denom);
p0 = 1/denom;
disp(['sum2=' num2str(sum2v) ' complnormcdf=' num2str(cmpl) ' p0=' num2str(p0) ' c=' num2str(y)]); %itc: HERE rm asap
end


% auxiliary functions


function y = sum1(r,R,T,t,l,sigma,m,IC,phat,p2,eps)
y = 0;
n=1;
last = 0; count = 0;
while true
    sum = quadvec(@(x) finner2(x,r,R,T,l,sigma,t,n,m,IC,phat,p2), 0, R-r);
    
    %itc20181006: the integral may run into instabilities resulting in zero
    %values because the shape of the function finner2 may come to look like
    %a gaussian with extrememly narrow spread, in which case, the numerical
    %integrator may come to think of it the result as zero; in such cases
    %we break up the interval [0,R-r] in many smaller intervals, integrate
    %over each one, and add up the results.
    if abs(sum)<10^-5 || ~isfinite(sum)
       % let's hope breaking up in smaller intervals will do the trick
       num_ints = 50;
       sz = (R-r)./num_ints;
       ll=0;
       ul=sz;
       for i=1:num_ints
           saux = quadvec(@(x) finner2(x,r,R,T,l,sigma,t,n,m,IC,phat,p2), ll, ul);
           %disp(['n=' num2str(n) ' ll=' num2str(ll) ' ul=' num2str(ul) ' saux=' num2str(saux) ' sum (so far)=' num2str(sum)]);
           if ~isnan(saux)
               sum = sum + saux;
           end
           ll=ul;
           ul = ul+sz;
       end
    end
    
    y = y+sum;
    last = last+sum;
    count = count + 1;
    if count == 4
        if abs(last/y) < eps || (abs(last)<10^-30 && abs(y)<10^-30)
            break;  % series converged
        else
            %disp(['sum1=' num2str(sum) ' y=' num2str(y) ' n=' num2str(n)]);
            count = 0;
            last = 0;
        end
    end
    n = n+1;
    if n==1000
        disp(['sum=' num2str(sum) ' last=' num2str(last) ' y=' num2str(y)]);
        error 'sum1: too many convolutions needed';
    end
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
%     disp(['sum2(r=' num2str(r) ',R=' num2str(R) ',T=' num2str(T) ',l=' num2str(l) ',s=' num2str(sigma) ',eps=' num2str(eps) ')']); %itc: HERE rm asap
%     disp(['calling integral(finner3(x,r,R,T,l,s,n=' num2str(n) ') from 0 to ' num2str(R-r)]); %itc: HERE rm asap
    sum = quadvec(@(x) finner3(x,r,R,T,l,sigma,n), 0, R-r);
%     disp(['sum=' num2str(sum)]); % itc: HERE rm asap
    %itc20181006: the integral may run into instabilities resulting in zero
    %values because the shape of the function finner3 may come to look like
    %a gaussian with extrememly narrow spread, in which case, the numerical
    %integrator may come to think of it the result as zero; in such cases
    %we break up the interval [0,R-r] in many smaller intervals, integrate
    %over each one, and add up the results.
    if abs(sum)<10^-5 || ~isfinite(sum)
       % let's hope breaking up in smaller intervals will do the trick
       num_ints = 50;
       sz = (R-r)./num_ints;
       ll=0;
       ul=sz;
       sum=0;
       for i=1:num_ints
           saux = quadvec(@(x) finner3(x,r,R,T,l,sigma,n), ll, ul);
           if ~isnan(saux)
            sum = sum + saux;
           end
           ll=ul;
           ul = ul+sz;
       end
    end
    
    y = y+sum;
    last = last+sum;
    count = count + 1;
    if count == 4
        if abs(last/y) < eps || (abs(last)<10^-30 && abs(y)<10^-30)
            break;  % series converged
        else
            %disp(['sum2=' num2str(sum) ' y=' num2str(y) ' n=' num2str(n)]);
            count = 0;
            last = 0;
        end
    end
    n = n+1;
    if n==1000
        disp(['sum2=' num2str(sum1) ' last=' num2str(last) ' y=' num2str(y)]);
        error 'sum2: too many convolutions needed';
    end
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
if isnan(xstn)
    xstn=0;
end
xstp = (x+m.*T)./(s.*sqrt(T));
if isnan(xstp)
    xstp=0;
end
t1 = (T-x./m-(s.^2)/(2.*m.^2)).*complnormcdf(xstn,0,1);
t2 = s.*sqrt(T)*normpdf(xstn,0,1)./m;
%t3 = (s.^2).*exp(2.*m.*x./(s.^2)).*complnormcdf(xstp,0,1)./(2.*m.^2);
%itc20181005: rearranged computation to avoid Inf*0 multiplication
cncdfp2 = (complnormcdf(xstp,0,1)./(2.*m.^2));
if abs(cncdfp2)<1.e-8 % avoid Inf*0 mult
    %t3=0;
    % compute exp(2mx/s^2)*complnormcdf(xstp) more accurately
    lnq1 = 2.*m.*x./(s.^2);
    lnq2 = log(cncdfp2);
    q = exp(lnq1+lnq2);
    t3 = (s.^2).*q./(2.*m.^2);
else
    t3 = cncdfp2.*exp(2.*m.*x./(s.^2)).*(s.^2);
end
y = t1+t2+t3;
end

function y=V1(x,T,m,s)
xstn = (x-m.*T)./(s.*sqrt(T));
if isnan(xstn)
    xstn=0;
end
xstp = (x+m.*T)./(s.*sqrt(T));
if isnan(xstp)
    xstp=0;
end
t1 = complnormcdf(xstn,0,1).*(T.^2 - (x./m).^2 - 2.*(s.^2).*x./(m.^3) - 3.*(s.^4)./(2.*m.^4))./2;
t2 = normpdf(xstn,0,1).*s.*sqrt(T).*(m.*T + 3.*s.*s./m + x)./(2.*m.^2);
%t3 = (s.^2).*(x - 3.*(s.^2)./(2.*m)).*exp(2.*m.*x./(s.^2)).*complnormcdf(xstp,0,1)./(2.*m.^3);
%itc20181005: rearranged computation to avoid Inf*0 multiplication
cncdfp2 = (complnormcdf(xstp,0,1)./(2.*m.^3));
if abs(cncdfp2)<1.e-8 % avoid Inf*0 mult
    %t3=0;
    lnq1 = 2.*m.*x./(s.^2);
    lnq2 = log(cncdfp2);
    q = exp(lnq1+lnq2);
    t3 = q.*(s.^2).*(x - 3.*(s.^2)./(2.*m));
else
    t3 = cncdfp2.*exp(2.*m.*x./(s.^2)).*(s.^2).*(x - 3.*(s.^2)./(2.*m));
end
y = t1 + t2 - t3;
end

function y=W1(x,T,m,s)
xstn = (x-m.*T)./(s.*sqrt(T));
if isnan(xstn)
    xstn=0;
end
xstp = (x+m.*T)./(s.*sqrt(T));
if isnan(xstp)
    xstp=0;
end
t1 = (s.^2).*(1 + m.*x./(s.^2)).*complnormcdf(xstn,0,1)./(m.^3);
t2 = 2.*s.*sqrt(T).*normpdf(xstn,0,1)./(m.^2);
%t3 = (x-(s.^2)./m).*exp(2.*m.*x./(s.^2)).*complnormcdf(xstp,0,1)./(m.^2);
%itc20181005: rearranged computation to avoid Inf*0 multiplication
cncdfm2 = (complnormcdf(xstp,0,1)./(m.^2));
if abs(cncdfm2)<1.e-8 % avoid Inf*0 mult
    %t3=0;
    lnq1 = 2.*m.*x./(s.^2);
    lnq2 = log(cncdfm2);
    q = exp(lnq1+lnq2);
    t3 = (x-(s.^2)./m).*q;
else
    t3=cncdfm2.*exp(2.*m.*x./(s.^2)).*(x-(s.^2)./m);
end
%disp(['W1: t1=' num2str(t1) ' t2=' num2str(t2) ' t3=' num2str(t3) ' (x=' num2str(x) ',T=' num2str(T) ',m=' num2str(m) ',s=' num2str(s) ')']);
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
% y = H(r+x,T,t,l,sigma,m,IC,phat,p2).*normpdf(R-r-x, n.*l.*T, sigma.*sqrt(n.*T));
% if ~isfinite(y)
%     disp(['y=' num2str(y) '=finner2(x=' num2str(x) ',r=' num2str(r) ',R=' num2str(R) ',T=' num2str(T) ',l=' num2str(l) ',sigma=' num2str(sigma) ',t=' num2str(t) ',n=' num2str(n) ',m=' num2str(m) ',IC=' num2str(IC) ',phat=' num2str(phat) ',p2=' num2str(p2) ')']);
%     error 'bad y';
% end
y = H(r+x,T,t,l,sigma,m,IC,phat,p2).*normpdf(R-r-x, n.*l.*T, sigma.*sqrt(n.*T));
% arg1=H(r+x,T,t,l,sigma,m,IC,phat,p2);
% log1 = log(arg1);
% arg2=normpdf(R-r-x, n.*l.*T, sigma.*sqrt(n.*T));
% log2 = log(arg2);
% %disp(['arg1=' num2str(arg1) ' log1=' num2str(log1) ' arg2=' num2str(arg2) ' log2=' num2str(log2)]);
% y = exp(log1+log2);
if ~isfinite(y)
    disp(['y=' num2str(y) '=finner2(x=' num2str(x) ',r=' num2str(r) ',R=' num2str(R) ',T=' num2str(T) ',l=' num2str(l) ',sigma=' num2str(sigma) ',t=' num2str(t) ',n=' num2str(n) ',m=' num2str(m) ',IC=' num2str(IC) ',phat=' num2str(phat) ',p2=' num2str(p2) ')']);
    error 'bad y';
end
end

function y = finner3(x,r,R,T,l,sigma,n)
y = n.*complnormcdf(x, l.*T, sigma.*sqrt(T)).*normpdf(R-r-x, (n-1).*l.*T, sigma.*sqrt((n-1).*T));
end

