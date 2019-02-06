function [sopt Qopt costs] = snQTCnormOptGraph_Q(Qmax,T,Kr,K0,L,mi,sigma,h,p, epsq)
% plot the min_{s}c(s,Q,T) cost of the (s,nQ,T) policy as a function of Q
% for given T, and compare with (R,T) policy

if nargin < 10
    epsq = 1.0;
end

% first figure out the cost of the (S,T) policy for the given T
%[sqt tunused c_ST] = STCnormOpt2(Kr,K0,L,mi,sigma,h,p,T);

Qlen = Qmax/epsq;

costs=1:Qlen;
Qi=1:Qlen;
hmax_Qi=1:Qlen;  % h(Q) -> the min_{s} val of the cost function H(s,Q,T;Kr+K0)
hmin_Qi=1:Qlen;  % h(T) -> as above but with Kr only
c_STi = 1:Qlen;
step_s = 0.001;
Qopt = -1;
copt=10.0^30;
s0=mi*(L+T);
for Q=1:Qlen
    Qi(Q)=(Q-1)*epsq + 1.e-6;  % 1.e-xxx is substitute for Q=0.
%     smin=-(mi+10.0*sigma)*(L+T);
%     smax=(mi+10.0*sigma)*(L+T);
%     sqt=lsqnonlin(@(x) aTLCeq(x,Qi(Q),T,L,mi,sigma,h,p), s0, smin, smax);
    sqt = findmins3(Qi(Q),T,L,mi,sigma,h,p,step_s,s0);
    % check soln
    v = aTLCeq(sqt,Qi(Q),T,L,mi,sigma,h,p);
    if abs(v)>1.e-3
        error(['findmins3() has error sqt=' num2str(sqt) ' v=' num2str(v)]);
    end
    s0 = sqt;
    % sanity check
    %if abs(aTLCeq(sqt,Qi(Q),T,L,mi,sigma,h,p)) > 1.e-5
    %    error(['for Q=' num2str(Qi(Q)) ' s=' num2str(sqt) ' is not optimized']);
    %end
    c = snQTCnorm(sqt,Qi(Q),T,Kr,K0,L,mi,sigma,h,p);
    hmax_Qi(Q) = snQTCnorm(sqt,Qi(Q),T,Kr+K0,0,L,mi,sigma,h,p);
    hmin_Qi(Q) = snQTCnorm(sqt,Qi(Q),T,Kr,0,L,mi,sigma,h,p);
    costs(Q)=c;
    if c < copt
       copt = c;
       sopt = sqt;
       Qopt=Qi(Q);
    end
end
hold on
plot(Qi(1:Qlen),hmin_Qi(1:Qlen),'g-');
hold off
hold on
plot(Qi(1:Qlen), hmax_Qi(1:Qlen), 'r-');
hold off
% for i=1:Qlen
%     c_STi(i)=c_ST;
% end
% hold on
% plot(Qi,c_STi,'b-');
% hold off
hold on
plot(Qi(1:Qlen),costs(1:Qlen),'k.-');
hold off
end

function sqt = findmins3(Q,T,L,mi,sigma,h,phat,eps_s,s0)
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
    [dir news cats] = detdir(s,Q,T,L,mi,sigma,h,phat,eps_s);
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

function [dir news cats] = detdir(s,Q,T,L,mi,sigma,h,phat,eps_s)
c = snQTCnorm(s,Q,T,0,0,L,mi,sigma,h,phat);
cats = c;
cup=snQTCnorm(s+eps_s,Q,T,0,0,L,mi,sigma,h,phat);
if c > cup
    dir=1;
    news=s;
    return;
else if c==cup
        s2=s;
        while snQTCnorm(s2,Q,T,0,0,L,mi,sigma,h,phat)==c
            s2=s2+eps_s;
        end
        cnew = snQTCnorm(s2,Q,T,0,0,L,mi,sigma,h,phat);
        if cnew < c
            dir=1;
            news=s2;
            return;
        else % cnew > c
            s2 = s;
            while snQTCnorm(s2,Q,T,0,0,L,mi,sigma,h,phat)==c
                s2 = s2-eps_s;
            end
            cnew = snQTCnorm(s2,Q,T,0,0,L,mi,sigma,h,phat);
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
cdown = snQTCnorm(s-eps_s,Q,T,0,0,L,mi,sigma,h,phat);
if c > cdown
    dir=-1;
    news=s;
    return;
else if c==cdown
        s2=s;
        while snQTCnorm(s2,Q,T,0,0,L,mi,sigma,h,phat)==c
            s2=s2-eps_s;
        end
        cnew = snQTCnorm(s2,Q,T,0,0,L,mi,sigma,h,phat);
        if cnew < c
            dir=-1;
            news=s2;
            return;
        else % cnew > c
            s2 = s;
            while snQTCnorm(s2,Q,T,0,0,L,mi,sigma,h,phat)==c
                s2 = s2+eps_s;
            end
            cnew = snQTCnorm(s2,Q,T,0,0,L,mi,sigma,h,phat);
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

