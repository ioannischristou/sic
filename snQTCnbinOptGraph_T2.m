function [sopt Qopt Topt copt] = snQTCnbinOptGraph_T2(Qmax,Tmax,Kr,K0,L,lamda,pl,h,p,epsq,epst,Qmin,Tmin,step_s)

mtdem = -lamda.*pl/((1-pl).*log(1-pl)); 

if nargin < 13
    Tmin = 5/mtdem;
end
if nargin < 12
    Qmin = 1;
end
if nargin < 11
    epst = 0.01;
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
copt = +Inf;
sqt = round(mtdem.*(L+Tmin));
for i=1:Tlen
%    hoptB=10.0^30;
%    hoptS=10.0^30;
    c_approx_Ti(i)=+Inf;
    cTi(i)=+Inf;
    T=Tmin+(i-1)*epst;
    Ti(i)=T;
    for j=1:Qlen
        Q=Qmin+(j-1)*epsq;
        Qi(j)=Q;
        s0=sqt;
        sqt = findmins3(Q,T,L,lamda,pl,h,p,step_s,s0);
        [c hm] = snQTCnbin2(sqt,Q,T,Kr,K0,L,lamda,pl,h,p);
        c_approx = snQTCnbinApprox2(sqt,Q,T,Kr,K0,L,lamda,pl,h,p);
        %disp(['for Q=' num2str(Q) ',T=' num2str(T) ': sqt=' num2str(sqt) ' c=' num2str(c) ' c_approx=' num2str(c_approx) ' hm=' num2str(hm) ' cur_copt=' num2str(copt)]);
        if c < cTi(i)
            disp(['T=' num2str(T) ' cTi=' num2str(cTi(i)) ' c=' num2str(c) ' @ Q=' num2str(Q) ',s=' num2str(sqt)]);
            %if 100*(cTi(i)-c)/c > 5.0 || (i>1 && Qopti(i-1)>=Q) % itc20181209: added tolerance check
                cTi(i)=c;
                %copt = c;
                Qopti(i)=Q;
                sopti(i)=sqt;
            %else
            %    disp(['ignoring bullshit difference...']);
            %end
        end
        if c_approx < c_approx_Ti(i)
            c_approx_Ti(i) = c_approx;
        end
        if c < copt
           sopt = sqt;
           Qopt = Q;
           Topt = T;
           copt = c;
        end
        % check if Q-search can stop
        B_L = hm;  % used to be hm + Kr./T;
        if B_L > cTi(i)  % itc20181211: used to be copt
                         % which would be correct for locating the global
                         % optimum, but incorrect for showing the graph
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
% finally, plot Qopti
% mTi = 1:Tlen;
% for i=1:Tlen
%     mTi(i)=mi.*Ti(i);
% end
figure()
hold on
plot(Ti,sopti,'g-');
hold off
hold on
plot(Ti,Qopti,'b');
hold off

end


function sqt = findmins3(Q,T,L,lamda,pl,h,phat,eps_s,s0)
%disp(['findmins3 Q=' num2str(Q) ' T=' num2str(T)]);
step = eps_s;
%itc20180910: niterbnd used to be 5
%niterbnd = 5;
niterbnd = 3;
mul = 2;
s = s0;
sqt = s;
cnt = niterbnd;
prevdir=0;
while true
    cnt = cnt-1;
    if cnt==0
        step=step*mul;
        %disp(['step=' num2str(step)]);
        cnt=niterbnd;
    end
    [dir news cats] = detdir(s,Q,T,L,lamda,pl,h,phat,eps_s);
    if news ~= s
        %disp(['s=' num2str(s) ' news=' num2str(news) ' cur_c=' num2str(cats)]);  % itc: HERE rm asap
    end
    if dir==0
        sqt = news;
        break;
    end
    if dir>0
        if prevdir<0 && step>eps_s
            step = step/mul;
            cnt=niterbnd;
            %disp(['step=' num2str(step)]);
        end
        s = s+step;  % itc: HERE this should probably go before if stmt
    else % dir<0
        s = s-step;
        if prevdir>0 && step>eps_s
            step = step/mul;
            cnt=niterbnd;
            %disp(['step=' num2str(step)]);
        end
    end
    %c2 = snQTCnorm2(s,Q,T,0,0,L,mi,sigma,h,phat); % itc: HERE rm asap 20181205
    %disp(['s=' num2str(s) ' dir=' num2str(dir) ' prevdir=' num2str(prevdir) ' cur_c=' num2str(c2)]);    
    prevdir=dir;
end
%disp(['findmins3() sqt=' num2str(sqt)]);
end

function [dir news cats] = detdir(s,Q,T,L,lamda,pl,h,phat,eps_s)
c = snQTCnbin(s,Q,T,0,0,L,lamda,pl,h,phat);
cats = c;
cup=snQTCnbin(s+eps_s,Q,T,0,0,L,lamda,pl,h,phat);
if c > cup
    dir=1;
    news=s;
    return;
else if c==cup
        s2=s;
        while snQTCnbin(s2,Q,T,0,0,L,lamda,pl,h,phat)==c
            s2=s2+eps_s;
        end
        cnew = snQTCnbin(s2,Q,T,0,0,L,lamda,pl,h,phat);
        if cnew < c
            dir=1;
            news=s2;
            return;
        else % cnew > c
            s2 = s;
            while snQTCnbin(s2,Q,T,0,0,L,lamda,pl,h,phat)==c
                s2 = s2-eps_s;
            end
            cnew = snQTCnbin(s2,Q,T,0,0,L,lamda,pl,h,phat);
            if cnew < c
                dir = -1;
                news=s2;
                return;
            else % cnew > c
                dir = 0;
                news = s;
                return;
            end
        end
    end
end
cdown = snQTCnbin(s-eps_s,Q,T,0,0,L,lamda,pl,h,phat);
if c > cdown
    dir=-1;
    news=s;
    return;
else if c==cdown
        s2=s;
        while snQTCnbin(s2,Q,T,0,0,L,lamda,pl,h,phat)==c
            s2=s2-eps_s;
        end
        cnew = snQTCnbin(s2,Q,T,0,0,L,lamda,pl,h,phat);
        if cnew < c
            dir=-1;
            news=s2;
            return;
        else % cnew > c
            s2 = s;
            while snQTCnbin(s2,Q,T,0,0,L,lamda,pl,h,phat)==c
                s2 = s2+eps_s;
            end
            cnew = snQTCnbin(s2,Q,T,0,0,L,lamda,pl,h,phat);
            if cnew < c
                dir = 1;
                news=s2;
                return;
            else % cnew > c
                dir = 0;
                news = s;
                return;
            end
        end
    end
end
% if we reach here we're optimal
dir=0;
news=s;
end

