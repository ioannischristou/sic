function [sopt Qopt Topt copt] = snQTCnormOptGraph_T3(Qmax,Tmax,Kr,K0,L,mi,sigma,h,p,epsq,epst,Qmin,Tmin,step_s)
if nargin < 14
    step_s = 1.e-6;
end
if nargin < 13
    Tmin = (3.5*sigma/mi)^2;
end
if nargin < 12
    Qmin = 1.e-8;
end
if nargin < 11
    epst = 0.1;
end
if nargin < 10
    epsq = 0.1;
end

p2=0;

Qlen = (Qmax-Qmin)/epsq;
Tlen = (Tmax-Tmin)/epst;

T_star_eoq = sqrt(2*K0.*(p+h)./(mi.*h.*p));
c_star_eoq = sqrt(2*K0.*mi.*p.*h./(p+h));

Qi=1:Qlen;
Ti=1:Tlen;
cTi=1:Tlen;
c_approx_Ti=1:Tlen;  % approximate cost using approximation of Porder
cSTi=1:Tlen;
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
        sqt = findmins3(Q,T,L,mi,sigma,h,p,step_s,s0,Qmin);
        if Q==Qmin
            cSTi(i)=snQTCnorm(sqt,Q,T,Kr+K0,0,L,mi,sigma,h,p,p2,Qmin);
        end
        c = snQTCnorm(sqt,Q,T,Kr,K0,L,mi,sigma,h,p,p2,Qmin);
        % check for errors
        cup = snQTCnorm(sqt+step_s,Q,T,Kr,K0,L,mi,sigma,p,p,p2,Qmin);
        if cup<c
            disp(['for Q=' num2str(Q) ',T=' num2str(T) ': sqt=' num2str(sqt) ' c=' num2str(c) ' c_approx=' num2str(c_approx)]);
            error(['cup=' num2str(cup)]);
        end
        cdown = snQTCnorm(sqt-step_s,Q,T,Kr,K0,L,mi,sigma,p,p,p2,Qmin);
        if cdown<c
            disp(['for Q=' num2str(Q) ',T=' num2str(T) ': sqt=' num2str(sqt) ' c=' num2str(c) ' c_approx=' num2str(c_approx)]);
            error(['cdown=' num2str(cdown)]);
        end        
        c_approx = snQTCnormApprox2(sqt,Q,T,Kr,K0,L,mi,sigma,h,p,p2,Qmin);
        disp(['for Q=' num2str(Q) ',T=' num2str(T) ': sqt=' num2str(sqt) ' c=' num2str(c) ' c_approx=' num2str(c_approx)]);
        if isnan(c) || isnan(c_approx)
            error('fucked up');
        end
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
        B_L = snQTCnorm(sqt,Q,T,Kr,0,L,mi,sigma,h,p,p2,Qmin);
        if B_L > cTi(i)  % used to be copt but then the plot would be incorrect
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
% figure()
% hold on
% plot(Ti,Qopti,'b');
% hold off
% hold on
% plot(Ti,sopti,'y');
% hold off
figure()
subplot(2,1,1);
hold on
plot(Ti,c_approx_Ti,'color',[0.8 0.8 0.8]);
hold off
hold on
plot(Ti,cTi, 'k.');
hold off
% hold on
% plot(Ti,cSTi,'b');
% hold off
hold on
plot(T_star_eoq,c_star_eoq,'r*');
hold off

% print out the difference between c and capprox
cdiff_Ti = 1:Tlen;
for i=1:Tlen
    cdiff_Ti(i)=cTi(i)-c_approx_Ti(i);
end
% figure()
subplot(2,1,2)
hold on
plot(Ti,cdiff_Ti,'r');
hold off

% figure()
% stdiffi=1:Tlen;
% for i=1:Tlen
%     stdiffi(i)=cTi(i)-cSTi(i);
% end
% hold on
% plot(Ti,stdiffi,'r.');
% hold off

% % finally, plot sopti and Qopti
% % mTi = 1:Tlen;
% % for i=1:Tlen
% %     mTi(i)=mi.*Ti(i);
% % end

figure()
% hold on
% plot(Ti,sopti,'g');
% hold off
hold on
plot(Ti,Qopti,'b');
hold off

end



function sqt = findmins3(Q,T,L,mi,sigma,h,phat,eps_s,s0,Qmin)
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
    [dir news cats] = detdir(s,Q,T,L,mi,sigma,h,phat,eps_s,Qmin);
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

function [dir news cats] = detdir(s,Q,T,L,mi,sigma,h,phat,eps_s,Qmin)
c = snQTCnorm(s,Q,T,0,0,L,mi,sigma,h,phat,0,Qmin);
cats = c;
cup=snQTCnorm(s+eps_s,Q,T,0,0,L,mi,sigma,h,phat,0,Qmin);
if c > cup
    dir=1;
    news=s;
    return;
else if c==cup
        s2=s;
        while snQTCnorm(s2,Q,T,0,0,L,mi,sigma,h,phat,0,Qmin)==c
            s2=s2+eps_s;
        end
        cnew = snQTCnorm(s2,Q,T,0,0,L,mi,sigma,h,phat,0,Qmin);
        if cnew < c
            dir=1;
            news=s2;
            return;
        else % cnew > c
            s2 = s;
            while snQTCnorm(s2,Q,T,0,0,L,mi,sigma,h,phat,0,Qmin)==c
                s2 = s2-eps_s;
            end
            cnew = snQTCnorm(s2,Q,T,0,0,L,mi,sigma,h,phat,0,Qmin);
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
cdown = snQTCnorm(s-eps_s,Q,T,0,0,L,mi,sigma,h,phat,0,Qmin);
if c > cdown
    dir=-1;
    news=s;
    return;
else if c==cdown
        s2=s;
        while snQTCnorm(s2,Q,T,0,0,L,mi,sigma,h,phat,0,Qmin)==c
            s2=s2-eps_s;
        end
        cnew = snQTCnorm(s2,Q,T,0,0,L,mi,sigma,h,phat,0,Qmin);
        if cnew < c
            dir=-1;
            news=s2;
            return;
        else % cnew > c
            s2 = s;
            while snQTCnorm(s2,Q,T,0,0,L,mi,sigma,h,phat,0,Qmin)==c
                s2 = s2+eps_s;
            end
            cnew = snQTCnorm(s2,Q,T,0,0,L,mi,sigma,h,phat,0,Qmin);
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

