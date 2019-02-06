function [Topt cTi] = snQTCnormGraph_TFixedsQ(s,Q,Tmax,Kr,K0,L,mi,sigma,h,p, epst)
if nargin < 11
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
sqt = s;
for i=1:Tlen
    Ti(i)=i*epst;
    T = Ti(i);
    if mi*sqrt(T)/sigma < 3
        thres = i;
    end
    c = snQTCnorm2(sqt,Q,T,Kr,K0,L,mi,sigma,h,p);
    cTi(i)=c;
    B0i(i)=mi*K0/Q + snQTCnorm2(sqt,Q,T,Kr,0,L,mi,sigma,h,p);
    B1i(i) = snQTCnorm2(sqt,Q,T,Kr+K0,0,L,mi,sigma,h,p);
    [st tt ct] = STCnormOpt2(Kr,K0,L,mi,sigma,h,p,T);
    B2i(i) = ct;
    if c < copt
        Topt = T;
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