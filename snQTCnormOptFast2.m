function [sopt Qopt Topt copt Qmax Tmax copt2] = snQTCnormOptFast2(Kr,K0,L,mi,sigma,h,p,epsq,epst,qnot,tnot,epss)
copt = 10.0^30;   % the optimal (s,nQ,T) value; 
copt2 = 10.0^30;  % the optimal (s,nQ,T) value s.t. Q>epsilon
stop_at_first_Q = false; 
if nargin <11
    tnot=((3.5*sigma)/mi)^2;
end
if nargin <10
    qnot=10^-8;
end
if nargin < 9  % time dimension precision
    epst = 0.1;
end
if nargin < 8  % batch size dimension precision
    epsq = 1;
end
if nargin < 12
    epss=0.01;
end
if nargin == 11 && tnot < 0
    stop_at_first_Q=true;
    tnot=((3*sigma)/mi)^2;
end
disp(['snQTCnormOptFast2: running with epsq=' num2str(epsq) ' epst=' num2str(epst) ' qnot=' num2str(qnot) ' tnot=' num2str(tnot) ' epss=' num2str(epss)]);
T=tnot;
Qmax=0;
%fudge = 0.0001;  % safety feature to account for issues with the derivation of P0

step_s = epss;

while T>0
    % Q=1e-3;
    % speed up T-search
    if T>=0.25 && epst<0.01
        epst = 0.01;
        disp(['T=' num2str(T) ' epst set to ' num2str(epst)]);
    else if T>=2 && epst<0.1
            epst = 0.1;
            disp(['T=' num2str(T) ' epst set to ' num2str(epst)]);
        else if T>=2.5 && epst<0.15
                epst=0.15;
                disp(['T=' num2str(T) ' epst set to ' num2str(epst)]);
            end
        end
    end
    Q=qnot;
    Hmdtprev = 10^25;
    Hm = 10^25;
    popprev = 10^25;
    s0 = mi.*(L+T);
    sopt_T = NaN;
    Qopt_T = NaN;
%    fq = 100;
    first_time = true;
    while Q>0
        if Q>2000
            error('Q got too big...');
        end
%         qfq = Q*fq;
%         diff = qfq - round(qfq);
%         if abs(diff) < 1.e-8  % itc20181212: speed up as Q gets bigger
%             prevepsq=epsq;
%             epsq = min(epsq*10,0.1);
%             if fq > 1
%                 fq = fq/10;
%             end
%             if prevepsq < 0.1
%                 disp(['epsq set to: ' num2str(epsq) ', Q=' num2str(Q)]);
%             end
%         end
        % let s' be the argmin. of c(s,Q,T)
%         s0=mi*(L+T);
        %smax=(mi+10.0*sigma)*(L+T);
        %smin = 0;
        %smax = (mi+50.0*sigma)*(L+T);
        %smin = -smax;
%         smax = 10^30;
%         smin = -10^30;
%         options = optimoptions(@lsqnonlin, 'FunctionTolerance', 1.e-8, 'TolX', 1.e-9, 'StepTolerance', 1.e-12);
%         [sqt val]=lsqnonlin(@(x) aTLCeq(x,Q,T,L,mi,sigma,h,p), s0, smin, smax, options);
%         if abs(val) > 10^-3
%             error(['h=' num2str(h) ' p=' num2str(p) ' Q=' num2str(Q) ' T=' num2str(T) ' s=' num2str(sqt) ' not solved val=' num2str(val)]);
%         end
        %sqt = findmins2(Q,T,L,mi,sigma,h,p,step_s,s0);
        sqt = findmins3(Q,T,L,mi,sigma,h,p,step_s,s0);
        s0 = sqt;
        %[c Hmdt] = snQTCnormM2(sqt,Q,T,Kr,K0,L,mi,sigma,h,p);
        [c Hmdt] = snQTCnorm(sqt,Q,T,Kr,K0,L,mi,sigma,h,p);
        %disp(['for Q=' num2str(Q) ',T=' num2str(T) ' sqt=' num2str(sqt) ' c=' num2str(c)]);  % itc: HERE rm asap
        if c < copt
           sopt = sqt;
           Qopt = Q;
           Topt = T;
           copt = c;
        end
        if Hmdt < Hm
            Hm = Hmdt;
            sopt_T = sqt;
            Qopt_T = Q;
        end
        if c < copt2 && Q>qnot
            copt2 = c;
        end
        % find the Hmdt = min_{s}H(s,Q,T;Kr)
        %sqt2 = lsqnonlin(@(x) aTLCeq(x,Q,T,L,mi,sigma,h,p), s0, smin, smax);
        %Hmdt = snQTCnorm(sqt,Q,T,Kr,0,L,mi,sigma,h,p);
        % b = mi*T+3*sigma*sqrt(T);
        %if Hmdt + K0*mi/b > copt && Hmdt > Hmdtprev
        if Hmdt > copt && Hmdt > Hmdtprev
            break;
        end
        % break if Po(Q,T) has passed its min, and -dH/dQ is less than
        % dPo/dQ
%         [pop hp passmin] = derivs(sqt,Q,T,K0,L,mi,sigma,h,p,popprev);
%         if passmin==1 && pop > -hp + fudge
%             break;
%         end
%         popprev = pop;
        Hmdtprev=Hmdt;
        if first_time
            Q=0;
            %disp(['T=' num2str(T) ' Q set to zero']);
            first_time=false;
        end
        Q=Q+epsq;
        if Qmax < Q
            Qmax=Q;
        end
        if stop_at_first_Q 
            break;
        end
    end
    % find the Hm = min_{s,Q}H(s,Q,T;Kr)
    %v0=[sqt2 Q];
    %vsol = fminsearch(@(v) aTLCeqH(v,T,Kr,L,mi,sigma,h,p), v0);
    %Hm = snQTCnorm(vsol(1), vsol(2), T, Kr,0,L,mi,sigma,h,p);
    % b = mi*T+3*sigma*sqrt(T);
    % if Hm + K0*mi/b > copt
    %Hm2 = snQTCnorm2(sqt,qnot,T,0,0,L,mi,sigma,h,p);
    if Hm > copt  % used to be Hm2 > copt
        break;
    end
    T=T+epst;
end
Tmax=T;
end


% function [Pop Hp passmin] = derivs(s,Q,T,K0,L,mi,sigma,h,p,popprev)
% hq=10^-7;
% l=mi;
% D=sigma^2;
% Pop = (K0/T)*(Porder(Q+hq,T,l,D)-Porder(Q-hq,T,l,D))/(2*hq);
% if Pop > popprev
%     passmin = 1;
% else
%     passmin = 0;
% end
% Hp = (snQTCnorm2(s,Q+hq,T,0,0,L,mi,sigma,h,p)-snQTCnorm2(s,Q-hq,T,0,0,L,mi,sigma,h,p))/(2*hq);
% end
% 
% 
% function y=Porder(Q,T,l,D)
% P01 = l*T*(normcdf((Q-l*T)/sqrt(D*T)))/Q;
% P02 = 1-normcdf((Q-l*T)/sqrt(D*T));
% P03 = normpdf((Q-l*T)/sqrt(D*T))*sqrt(D*T)/Q;
% y=P01+P02-P03;
% end
% 
% 
% function sqt = findmins2(Q,T,L,mi,sigma,h,phat,step_s,s0)
% s = s0;
% sqt = s;
% c = snQTCnorm2(s,Q,T,0,0,L,mi,sigma,h,phat);
% cs = snQTCnorm2(s+step_s,Q,T,0,0,L,mi,sigma,h,phat);
% if c < cs 
%    % keep decreasing s
%    while true
%        s = s-step_s;
%        cs = snQTCnorm2(s,Q,T,0,0,L,mi,sigma,h,phat);
%        if cs >= c
%            sqt = s+step_s;
%            break;
%        end
%        c = cs;
%        sqt = s;
%    end
% else if c > cs
%         % keep increasing s
%         while true
%             s = s+step_s;
%             cs = snQTCnorm2(s,Q,T,0,0,L,mi,sigma,h,phat);
%             if cs >= c
%                 sqt = s-step_s;
%                 break;
%             end
%             c = cs;
%             sqt = s;
%         end
%     end
% end
% end


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

