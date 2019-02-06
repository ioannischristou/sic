function [sopt Qopt Topt costs] = snQTCpoissonOptGraph_T2(Tmax,Kr,K0,L,lamda,h,p, epsq, epst, tnot, p2, estop)
if nargin < 12
    estop=0;
end
if nargin < 11
    p2 = 0;
end
if nargin < 8
    epsq = 1;
end
if nargin < 9
    epst = 1;
end
if nargin < 10
    tnot = epst;
end
%Qlen = Qmax/epsq;
Tlen = (Tmax-tnot)/epst;
costs=1:Tlen;
%Qi=1:Qlen;
Ti=1:Tlen;
cTi=1:Tlen;
Qopti=1:Tlen;
sopti=1:Tlen;
hmax_Ti=1:Tlen;  % h(T) -> the min_{s,Q} val of the cost function H(s,Q,T;Kr+K0)
hmin_Ti=1:Tlen;  % h(T) -> as above but with Kr only 
copt_tot = 10.0^30;
mi = lamda;
%sigma = sqrt(lamda);
s0 = 1;  % initial estimate of optimal reorder point s
for i=1:Tlen
    copt=10.0^30;
    hoptB=10.0^30;
    hoptS=10.0^30;
    Ti(i)=(i-1)*epst+tnot;
    T = Ti(i);
    disp(['T=' num2str(T)]);
    j=1;
    while true
        j = j+1;
        Q = (j-1)*epsq;
        %Qi(j)=Q;
        %s0=mi*(L+T);
        %smin=-ceil(s0);
        % compute optimal s
        %sqt = findmins(Q,T,L,lamda,h,p,p2,smin);
        sqt = findmins2(Q,T,L,lamda,h,p,p2,s0);
        s0 = sqt;  % update s0 estimate for next Q
        c = snQTCpoisson(sqt,Q,T,Kr,K0,L,lamda,h,p,p2);
        hbig = snQTCpoisson(sqt,Q,T,Kr+K0,0,L,lamda,h,p,p2);
        hsmall = snQTCpoisson(sqt,Q,T,Kr,0,L,lamda,h,p,p2);
        %costs(j,i)=c;
        if c < copt
           cTi(i)=c;
           copt = c;
           Qopti(i)=Q;
           sopti(i)=sqt;
           costs(i)=c;
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
           disp(['found better value ' num2str(copt_tot) ' at Q=' num2str(Qopt) ' T=' num2str(Topt)]);
        end
        if hsmall > copt
            break;
        end
        if estop==1
            break;  % run the (S,T) policy Q=1
        end
    end
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
plot(Ti(1:Tlen), hmax_Ti(1:Tlen), 'r-');
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


function [sqt c] = findmins(Q,T,L,lamda,h,phat,p,smin)
sqt = smin;
c = snQTCpoisson(sqt,Q,T,0,0,L,lamda,h,phat,p);
s = sqt;
while true 
    s = s+1;
    cs = snQTCpoisson(s,Q,T,0,0,L,lamda,h,phat,p);
    if cs>=c
        break;
    end
    c = cs;
    sqt = s;
end
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