function y = sSTCpoissonGraph_S(s,Smax,T,Kr,K0,L,lamda,h,p,eS,Snot)
if nargin < 10
    eS = 1;
end
if nargin < 11
    Snot = s;
end
if Snot < s
    Snot = s;
end

Slen = (Smax-Snot)/eS;
Si = 1:Slen;
ci = 1:Slen;

Si(1)=Snot;
ci(1)=sSTCpoisson(s,Si(1),T,Kr,K0,L,lamda,h,p);
for i=2:Slen
    Si(i) = Si(i-1) + eS;
    ci(i) = sSTCpoisson(s,Si(i),T,Kr,K0,L,lamda,h,p);
    disp(['c(' num2str(Si(i)) ')=' num2str(ci(i))]);
end

hold on
plot(Si,ci);
hold off
y = ci;
end