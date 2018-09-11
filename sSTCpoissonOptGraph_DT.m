function [Di Ti costs] = sSTCpoissonOptGraph_DT(Dmin,Dmax,Tmin,Tmax,Kr,K0,L,lamda,h,p,epsd,epst)
if nargin < 11
    epsd = 1;
end
if nargin < 12
    epst = 1;
end
Di = Dmin:epsd:Dmax;
Ti = Tmin:epst:Tmax;
Dlen = numel(Di);
Tlen = numel(Ti);
costs = ones(Dlen,Tlen);
s0 = 1;  % init guess
for i=1:Tlen
    T = Ti(i);
    for j=1:Dlen
        D = Di(j);
        [s S c] = sSTCpoissonOpt_FixedDT(D,T,Kr,K0,L,lamda,h,p,0,s0);
        s0 = S;
        costs(j,i)=c;
    end
end
surf(Ti,Di,costs);
colormap hsv;
end