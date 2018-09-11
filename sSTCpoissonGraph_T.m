function y=sSTCpoissonGraph_T(Tlen,si,Si,Kr,K0,L,lamda,h,p,p2,epst,tnot)
if nargin < 10
    p2=0;
end
if nargin < 11
    epst = 1;
end
if nargin < 12
    tnot = 10^-6;
end
ci = 1:Tlen;
Ti = 1:Tlen;
for i=1:Tlen
    Ti(i) = tnot + (i-1)*epst;
    ci(i)=sSTCpoisson(si(i),Si(i),Ti(i),Kr,K0,L,lamda,h,p,p2);
end

plot(Ti,ci);
y=ci;
end