function y=finner2plot(r,R,T,l,sigma,t,n,m,IC,phat,p2,stepsize)
xi=0:stepsize:(R-r);
yi=1:numel(xi);
zi=1:numel(xi);
for i=1:numel(xi)
    x = xi(i);
    yi(i)=finner2(x,r,R,T,l,sigma,t,n,m,IC,phat,p2);
    if x>0
        zi(i)=quadvec(@(u) finner2(u,r,R,T,l,sigma,t,n,m,IC,phat,p2), 0, x);
    else
        zi(i)=0;
    end
end
hold on
plot(xi,yi,'k-');
hold off
hold on
plot(xi,zi,'b-');
end


function y=finner2(x,r,R,T,l,sigma,t,n,m,IC,phat,p2)
%y = H(r+x,T,t,l,sigma,m,IC,phat,p2).*normpdf(R-r-x, n.*l.*T, sigma.*sqrt(n.*T));
arg1=H(r+x,T,t,l,sigma,m,IC,phat,p2);
log1 = log(arg1);
arg2=normpdf(R-r-x, n.*l.*T, sigma.*sqrt(n.*T));
log2 = log(arg2);
disp(['arg1=' num2str(arg1) ' log1=' num2str(log1) ' arg2=' num2str(arg2) ' log2=' num2str(log2)]);
%ok, it's the normpdf that comes up zero, so print out args
%apdfx=R-r-x;
%apdfm=n.*l.*T;
%apdfs=sigma.*sqrt(n.*T);
%disp(['argnormpdf_x=' num2str(apdfx) ' argnormpdf_m=' num2str(apdfm) ' argnormpdf_s=' num2str(apdfs)]);
y = exp(log1+log2);
if ~isfinite(y)
    disp(['y=' num2str(y) '=finner2(x=' num2str(x) ',r=' num2str(r) ',R=' num2str(R) ',T=' num2str(T) ',l=' num2str(l) ',sigma=' num2str(sigma) ',t=' num2str(t) ',n=' num2str(n) ',m=' num2str(m) ',IC=' num2str(IC) ',phat=' num2str(phat) ',p2=' num2str(p2) ')']);
    error 'bad y';
end
end


function y = H(rpj,T,t,l,sigma,m,IC,phat,p2)
z = 0;
if nargin == 9 && p2~=0
    z = p2.*eP(rpj,T,l,sigma,t);
end
y = IC.*T.*(rpj - m - l.*T./2) + (IC+phat).*bP(rpj,T,l,sigma,t) + z;
%disp(['H(' num2str(rpj) ')=' num2str(y)]);
end


function y = complnormcdf(x,m,s)
y=1.-normcdf(x,m,s);
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
if abs(cncdfp2)<1.e-24 % avoid Inf*0 mult
    %t3=0;
    % compute exp(2mx/s^2)*complnormcdf(xstp) more accurately
    lnq1 = 2.*m.*x./(s.^2);
    lnq2 = log(cncdfp2);
    q = exp(lnq1+lnq2);
    if ~isfinite(q) 
        disp(['lnq1=' num2str(lnq1) ' lnq2= ' num2str(lnq2)]);
        error 'non-finite q';
    end
    t3 = (s.^2).*q./(2.*m.^2);
else
    t3 = cncdfp2.*exp(2.*m.*x./(s.^2)).*(s.^2);
    if ~isfinite(t3)
        disp(['t3=' num2str(t3)]);
        error 'non-finite t3';
    end
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
if abs(cncdfp2)<1.e-24 % avoid Inf*0 mult
    %t3=0;
    lnq1 = 2.*m.*x./(s.^2);
    lnq2 = log(cncdfp2);
    q = exp(lnq1+lnq2);
    t3 = q.*(s.^2).*(x - 3.*(s.^2)./(2.*m));
    if ~isfinite(t3)
        disp(['lnq1=' num2str(lnq1) ' lnq2=' num2str(lnq2)]);
        error 'non-finite t3';
    end
else
    t3 = cncdfp2.*exp(2.*m.*x./(s.^2)).*(s.^2).*(x - 3.*(s.^2)./(2.*m));
end
y = t1 + t2 - t3;
end

function y=W1(x,T,m,s)
disp(['calling W1(x=' num2str(x) ',T=' num2str(T) ',m=' num2str(m) ',s=' num2str(s) ')']);
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
if abs(cncdfm2)<1.e-24 % avoid Inf*0 mult
    %t3=0;
    lnq1 = 2.*m.*x./(s.^2);
    lnq2 = log(cncdfm2);
    q = exp(lnq1+lnq2);
    t3 = (x-(s.^2)./m).*q;
    if ~isfinite(t3)
       disp(['lnq1=' num2str(lnq1) ' lnq2=' num2str(lnq2)]);
       error 'non-finite t3';
    end
else
    t3=cncdfm2.*exp(2.*m.*x./(s.^2)).*(x-(s.^2)./m);
    if ~isfinite(t3)
        exp2_term = 2.*m.*x./(s.^2);
        exp2 = exp(exp2_term);
        last = (x-(s.^2)./m);
       disp(['t3=' num2str(t3) ' xstp=' num2str(xstp) ' cncdfm2=' num2str(cncdfm2) ' exp2_term=' num2str(exp2_term) ' exp2=' num2str(exp2) ' last=' num2str(last)]);
       error 'non-finite t3';
    end
end
%disp(['W1: t1=' num2str(t1) ' t2=' num2str(t2) ' t3=' num2str(t3) ' (x=' num2str(x) ',T=' num2str(T) ',m=' num2str(m) ',s=' num2str(s) ')']);
y = t1 - t2 + t3;
end

