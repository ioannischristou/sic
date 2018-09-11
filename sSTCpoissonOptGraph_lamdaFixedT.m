function y = sSTCpoissonOptGraph_lamdaFixedT(T,lmin,lmax,Kr,K0,L,h,p,p2,epsl)
if nargin<10
    epsl=1;
end
li=lmin:epsl:lmax;
li_len=numel(li);
ci = 1:li_len;
for i=1:li_len
    [s S c] = sSTCpoissonOpt_FixedT(T,Kr,K0,L,li(i), h,p,p2);
    ci(i)=c;
end
hold on
plot(li,ci,'b');
hold off
y=ci;
end

