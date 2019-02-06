function [sopt Qopt Topt costs] = snQTCpoissonOptGraph_T(Qmax,Tmax,Kr,K0,L,lamda,h,p, epsq, epst)
if nargin < 9
    epsq = 1;
end
if nargin < 10
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
mi = lamda;
sigma = sqrt(lamda);
s0=1;
s1=1;
for i=1:Tlen
    copt=10.0^30;
    hoptB=10.0^30;
    hoptS=10.0^30;
    Ti(i)=(i-1)*epst+1.e-6;  % compute the <nQ,r> system as well for i=1
    T = Ti(i);
    for j=1:Qlen
        Q = (j-1)*epsq + 1;
        Qi(j)=Q;
        %s0=mi*(L+T);
        %smin=0.0;
        %smax=(mi+10.0*sigma)*(L+T);
        %[sqt c exitflag] = fmincon(@(x) snQTCpoisson(x,Q,T,Kr,K0,L,lamda,h,p), s0, [],[],[],[], smin, smax);
        %if exitflag<=0
        %    error('optimization in snQTpoissonOptFast() failed');
        %end        
        sqt = findmins2(Q,T,L,lamda,h,p,0,s0);  % 0 is the p2 value
        if Q==1
            s1 = sqt;
        end
        c = snQTCpoisson(sqt,Q,T,Kr,K0,L,lamda,h,p);
        s0 = sqt;
        hbig = snQTCpoisson(sqt,Q,T,Kr+K0,0,L,lamda,h,p);
        hsmall = snQTCpoisson(sqt,Q,T,Kr,0,L,lamda,h,p);
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
    disp(['T=' num2str(T)]);
end
% re-size hmax_Ti
%hmax_Ti(1)=hmax_Ti(3);
%hmax_Ti(2)=hmax_Ti(3);
% end re-sizing
%cost = copt;
hold on
plot(Ti,hmin_Ti,'g-');
hold off
hold on
plot(Ti(5:Tlen), hmax_Ti(5:Tlen), 'r-');
hold off
hold on
plot(Ti,Qopti,'b-');
hold off
%hold on
%plot(Ti(1:Tlen),sopti(1:Tlen),'color',[0.1 0.1 0.1]);
%hold off
hold on
plot(Ti,cTi,'k.-');
hold off
end

function sqt = findmins2(Q,T,L,lamda,h,phat,p,s0)
s = s0;
sqt = s;
c = snQTCpoisson(s,Q,T,0,0,L,lamda,h,phat,p);
cs = snQTCpoisson(s+1,Q,T,0,0,L,lamda,h,phat,p);
if c < cs 
   % keep decreasing s
   while true
       s = s-1;
       cs = snQTCpoisson(s,Q,T,0,0,L,lamda,h,phat,p);
       if cs >= c
           sqt = s+1;
           break;
       end
       c = cs;
       sqt = s;
   end
else if c > cs
        % keep increasing s
        while true
            s = s+1;
            cs = snQTCpoisson(s,Q,T,0,0,L,lamda,h,phat,p);
            if cs >= c
                sqt = s-1;
                break;
            end
            c = cs;
            sqt = s;
        end
    end
end
end

