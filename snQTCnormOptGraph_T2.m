function [sopt Qopt Topt copt] = snQTCnormOptGraph_T2(Qmax,Tmax,Kr,K0,L,mi,sigma,h,p,epsq,epst,Qmin,Tmin,step_s)
if nargin < 13
    Tmin = (3.5*sigma/mi)^2;
end
if nargin < 12
    Qmin = 1.e-6;
end
if nargin < 11
    epst = 0.1;
end
if nargin < 10
    epsq = 1;
end
% last default is necessarily done here
if nargin < 14
    step_s = epsq;
end

Qlen = (Qmax-Qmin)/epsq;
Tlen = (Tmax-Tmin)/epst;

Qi=1:Qlen;
Ti=1:Tlen;
cTi=1:Tlen;
c_approx_Ti=1:Tlen;  % approximate cost using approximation of Porder
Qopti=1:Tlen;
sopti=1:Tlen;
%hmax_Ti=1:Tlen;  % h(T) -> the min_{s,Q} val of the cost function H(s,Q,T;Kr+K0)
%hmin_Ti=1:Tlen;  % h(T) -> as above but with Kr only 
copt = 10.0^30;
sqt = mi.*(L+Tmin);
for i=1:Tlen
%    hoptB=10.0^30;
%    hoptS=10.0^30;
    c_approx_Ti(i)=10.0^30;
    cTi(i)=10.0^30;
    T=Tmin+(i-1)*epst;
    Ti(i)=T;
    for j=1:Qlen
        Q=Qmin+(j-1)*epsq;
        Qi(j)=Q;
        s0=sqt;
        %smin=1.0;
        %smin=-(mi+100*sigma).*(L+T);
        %smax=(mi+10.0*sigma)*(L+T);
        %smax=(mi+100.0*sigma).*(L+T);
        %sqt=lsqnonlin(@(x) aTLCeq(x,Q,T,L,mi,sigma,h,p), s0, smin, smax);
        sqt = findmins3(Q,T,L,mi,sigma,h,p,step_s,s0);
        c = snQTCnorm(sqt,Q,T,Kr,K0,L,mi,sigma,h,p);
        c_approx = snQTCnormApprox2(sqt,Q,T,Kr,K0,L,mi,sigma,h,p);
        disp(['for Q=' num2str(Q) ',T=' num2str(T) ': sqt=' num2str(sqt) ' c=' num2str(c) ' c_approx=' num2str(c_approx)]);
        %hbig = snQTCnorm2(sqt,Q,T,Kr+K0,0,L,mi,sigma,h,p);
        %hsmall = snQTCnorm2(sqt,Q,T,Kr,0,L,mi,sigma,h,p);
        if c < cTi(i)
           cTi(i)=c;
           %copt = c;
           Qopti(i)=Q;
           sopti(i)=sqt;
        end
        if c_approx < c_approx_Ti(i)
            c_approx_Ti(i) = c_approx;
        end
%         if hbig < hoptB
%             hmax_Ti(i)=hbig;
%             hoptB = hbig;
%         end
%         if hsmall < hoptS
%             hmin_Ti(i) = hsmall;
%             hoptS = hsmall;
%         end
        if c < copt
           sopt = sqt;
           Qopt = Q;
           Topt = T;
           copt = c;
        end
        % check if Q-search can stop
        B_L = snQTCnorm(sqt,Q,T,Kr,0,L,mi,sigma,h,p);
        if B_L > cTi(i)
            break;
        end
    end
end
% re-size hmax_Ti
%hmax_Ti(1)=hmax_Ti(3);
%hmax_Ti(2)=hmax_Ti(3);
% end re-sizing
%cost = copt;
%hold on
%plot(Ti,cTi,'k.-');
%hold off
% hold on
% plot(Ti,hmin_Ti,'g*-');
% hold off
% hold on
% plot(Ti(2:Tmax), hmax_Ti(2:Tmax), 'r*-');
% hold off
figure()
% hold on
% plot(Ti,Qopti,'b');
% hold off
% hold on
% plot(Ti,sopti,'y');
% hold on
plot(Ti,c_approx_Ti,'color',[0.8 0.8 0.8]);
hold off
hold on
plot(Ti,cTi, 'k.');
hold off
% print out the difference between c and capprox
cdiff_Ti = 1:Tlen;
for i=1:Tlen
    cdiff_Ti(i)=cTi(i)-c_approx_Ti(i);
end
figure()
hold on
plot(Ti,cdiff_Ti,'r');
hold off
% finally, plot mi*Ti and Qopti
mTi = 1:Tlen;
for i=1:Tlen
    mTi(i)=mi.*Ti(i);
end
figure()
hold on
plot(Ti,mTi,'k.');
hold off
hold on
plot(Ti,Qopti,'b');
hold off

end


function sqt = findmins2(Q,T,L,mi,sigma,h,phat,step_s,s0)
s = s0;
sqt = s;
c = snQTCnorm2(s,Q,T,0,0,L,mi,sigma,h,phat);
cs = snQTCnorm2(s+step_s,Q,T,0,0,L,mi,sigma,h,phat);
if c < cs 
   % keep decreasing s
   while true
       s = s-step_s;
       cs = snQTCnorm2(s,Q,T,0,0,L,mi,sigma,h,phat);
       if cs >= c
           sqt = s+step_s;
           break;
       end
       c = cs;
       sqt = s;
   end
else if c > cs
        % keep increasing s
        while true
            s = s+step_s;
            cs = snQTCnorm2(s,Q,T,0,0,L,mi,sigma,h,phat);
            if cs >= c
                sqt = s-step_s;
                break;
            end
            c = cs;
            sqt = s;
        end
    end
end
end

