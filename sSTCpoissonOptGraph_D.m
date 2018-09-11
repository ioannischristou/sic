function y=sSTCpoissonOptGraph_D(Dmax,T,Kr,K0,L,lamda,h,p,p2,epsd,dnot)
if nargin < 9
    p2 = 0;
end
if nargin < 10
    epsd = 1;
end
if nargin < 11
    dnot = 1;
end

Dlen = (Dmax-dnot)/epsd;

ci = 1:Dlen;
li = 1:Dlen;
ui = 1:Dlen;
Di = 1:Dlen;
S0 = 1;
S2 = 1;
S3 = 1;
for i=1:Dlen
    D = dnot + (i-1)*epsd;
    Di(i)=D;
    [s S c]=sSTCpoissonOpt_FixedDT(D,T,Kr,K0,L,lamda,h,p,p2, S0);
    if i==1
        S2 = S;
        S3 = S;
    end
    [s2 S2n l]=sSTCpoissonOpt_FixedDT(D,T,Kr,0,L,lamda,h,p,p2, S2);
    [s3 S3n u]=sSTCpoissonOpt_FixedDT(D,T,Kr+K0,0,L,lamda,h,p,p2, S3);
    ci(i)=c;
    li(i)=l;
    ui(i)=u;
    S0 = S;
    S2 = S2n;
    S3 = S3n;
    disp(['c(' num2str(D) ')=' num2str(c)]);
end

hold on
plot(Di,ci,'b');
plot(Di,li,'g');
plot(Di,ui,'r');
hold off
y = ci;

end


