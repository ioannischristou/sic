function [sopt Topt cTi] = snQTCnormOptGraph_TFixedQ(Q,Tmax,Kr,K0,L,mi,sigma,h,p, epst)
if nargin < 10
    epst = 1;
end
Tlen = Tmax/epst;
Ti=1:Tlen;
cTi=1:Tlen;
B0i = 1:Tlen;
B1i = 1:Tlen;
B2i = 1:Tlen;
copt = 10.0^30;
thres = 1;
for i=1:Tlen
    Ti(i)=i*epst;
    T = Ti(i);
    if mi*sqrt(T)/sigma < 3
        thres = i;
    end
    s0=mi*(L+T);
    smin=0.0;
    smax=(mi+10.0*sigma)*(L+T);
    sqt=lsqnonlin(@(x) aTLCeq(x,Q,T,L,mi,sigma,h,p), s0, smin, smax);
    c = snQTCnorm(sqt,Q,T,Kr,K0,L,mi,sigma,h,p);
    cTi(i)=c;
    B0i(i)=mi*K0/Q + snQTCnorm(sqt,Q,T,Kr,0,L,mi,sigma,h,p);
    B1i(i) = snQTCnorm(sqt,Q,T,Kr+K0,0,L,mi,sigma,h,p);
    [st tt ct] = STCnormOpt2(Kr,K0,L,mi,sigma,h,p,T);
    B2i(i) = ct;
    if c < copt
        Topt = T;
        sopt = sqt;
        copt = c;
    end
end
hold on
plot(Ti(thres:Tlen),B0i(thres:Tlen),'g-');
hold off
hold on
plot(Ti(thres:Tlen),B1i(thres:Tlen),'r-');
hold off
hold on
plot(Ti(thres:Tlen),B2i(thres:Tlen),'b-');
hold off
hold on
plot(Ti(thres:Tlen),cTi(thres:Tlen),'k-');
hold off

end