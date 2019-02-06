function [sopt Qopt Topt costs] = snQTCnorm2OptGraph_T(Qmax,Tmax,Kr,K0,L,mi,sigma,h,p, epsq, epst)
if nargin < 10
    epsq = 1;
end
if nargin < 11
    epst = 1;
end
Qlen = Qmax/epsq;
Tlen = Tmax/epst;
costs=ones(Qlen,Tlen);
Qi=1:Qlen;
Ti=1:Tlen;
cTi=1:Tlen;
Qopti=1:Tlen;
sopti=1:Tlen;
hmax_Ti=1:Tlen;  % h(T) -> the min_{s,Q} val of the cost function H(s,Q,T;Kr+K0)
hmin_Ti=1:Tlen;  % h(T) -> as above but with Kr only 
copt_tot = 10.0^30;
for i=1:Tlen
    copt=10.0^30;
    hoptB=10.0^30;
    hoptS=10.0^30;
    Ti(i)=i*epst;
    T = Ti(i);
    for j=1:Qlen
        Q = (j-1)*epsq + 1.e-6;
        Qi(j)=Q;
        s0=mi*(L+T);
        %smin=1.0;
        %smax=(mi+20.0*sigma)*(L+T);
        smax=10^5;
        smin=-smax;
        sqt=lsqnonlin(@(x) aTLCeq(x,Q,T,L,mi,sigma,h,p), s0, smin, smax);
        val = aTLCeq(sqt,Q,T,L,mi,sigma,h,p);
        if abs(val)>10^-3
            error('not solved');
        end
        c = snQTCnorm2(sqt,Q,T,Kr,K0,L,mi,sigma,h,p);
        hbig = snQTCnorm2(sqt,Q,T,Kr+K0,0,L,mi,sigma,h,p);
        hsmall = snQTCnorm2(sqt,Q,T,Kr,0,L,mi,sigma,h,p);
        costs(j,i)=c;
        if c < copt
           cTi(i)=c;
           copt = c;
           Qopti(i)=Q;
           sopti(i)=sqt;
        end
        if hbig < hoptB
            hmax_Ti(i)=hbig;
            hoptB = hbig;
        end
        if hsmall < hoptS
            hmin_Ti(i) = hsmall;
            hoptS = hsmall;
        end
        if c < copt_tot
           sopt = sqt;
           Qopt = Q;
           Topt = T;
           copt_tot = c;
        end
    end
end
% re-size hmax_Ti
%hmax_Ti(1)=hmax_Ti(3);
%hmax_Ti(2)=hmax_Ti(3);
% end re-sizing
%cost = copt;
% hold on
% plot(Ti,hmin_Ti,'g-');
% hold off
% hold on
% plot(Ti(6:Tlen), hmax_Ti(6:Tlen), 'r-');
% hold off
hold on
plot(Ti,Qopti,'b-');
hold off
%hold on
%plot(Ti(1:Tlen),sopti(1:Tlen),'color',[0.1 0.1 0.1]);
%hold off
hold on
plot(Ti,cTi,'g.-');
hold off
end