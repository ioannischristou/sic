function [sopt qopt topt copt] = snQTCnormOptFastApprox3(Kr, K0, L, mi, sigma, h, p, epsq, epst, qnot, tnot, epss)
% the function computes an approximate optimum point of the (r,nQ,T) policy
% when demands are normal, according to Lagodimos. 
if nargin < 11
    tnot = ((3.5*sigma)./mi).^2;
end
if nargin < 10
    qnot = 10^-8;
end
if nargin < 9
    epst = 0.1;
end
if nargin < 8
    epsq = 1.0;
end
if nargin < 12
    epss = 0.01;
end

%initialize return variables
sopt=NaN;
qopt=NaN;
topt=NaN;
copt=+Inf;

% compute EOQ quantities
T_star_eoq = sqrt(2*K0.*(p+h)./(mi.*h.*p));
Q_star_eoq = sqrt(2*K0.*mi.*(p+h)./(h.*p));

step_s = epss;

% compute optimal (S,T) 
start_T = max(ceil(T_star_eoq/epst)*epst,tnot);  % used to be T_star_eoq
[s_opt_st, t_opt_st, c_opt_st] = STCnormOptFastApprox(start_T,Kr,K0,L,mi,sigma,h,p,epst,qnot);


Tmax = 100;  % max limit for T-search
Qmax = 2000; % max limit for Q-search

T=tnot;
T_limit_st = c_opt_st-Kr./t_opt_st;

% outer-loop: T-search
while T<Tmax
    % short-cut heuristic when Kr=0
    if Kr==0 && T>tnot
        break;
    end
    % speed up T-search
    if T>=0.25 && epst<0.01
        epst = 0.01;
    else if T>=2 && epst<0.1
            epst = 0.1;
        else if T>=2.5 && epst<0.15
                epst=0.15;
            end
        end
    end
    Q=ceil(Q_star_eoq/epsq)*epsq;  % used to be Q_star_eoq
    % inner-loop: Q-search
    s0=mi*(L+T);    
    c_prev = +Inf;
%    fq = 100;
    while Q<Qmax
%         qfq = Q*fq;
%         diff = abs(qfq - round(qfq));
%         if diff < 1.e-8  % itc20181212: speed up as Q gets bigger
%             prevepsq = epsq;
%             epsq = min(epsq*10,0.1);
%             if fq > 1
%                 fq = fq/10;
%             end
%             if prevepsq < 0.1
%                 disp(['epsq set to: ' num2str(epsq) ', Q=' num2str(Q)]);
%             end
%         end
        %smax=(mi+10.0*sigma)*(L+T);
        %smin = 0;
        %smax = (mi+50.0*sigma)*(L+T);
        %smin = -smax;
%         smax = 10^30;
%         smin = -10^30;
%         [sqt, val]=lsqnonlin(@(x) aTLCeq(x,Q,T,L,mi,sigma,h,p), s0, smin, smax);
%         if abs(val) > 10^-4
%             error(['h=' num2str(h) ' p=' num2str(p) ' Q=' num2str(Q) ' T=' num2str(T) ' s=' num2str(sqt) ' not solved val=' num2str(val)]);
%         end
        %sqt = findmins2(Q,T,L,mi,sigma,h,p,step_s,s0);
        sqt = findmins3(Q,T,L,mi,sigma,h,p,step_s,s0);
        % sqt is the current cost minimizer given (Q,T)
        c = snQTCnormApprox2(sqt,Q,T,Kr,K0,L,mi,sigma,h,p);
        if c<copt
            sopt=sqt;
            qopt=Q;
            topt=T;
            copt=c;
        end
        s0=sqt;
        if c_prev < c
            break;
        end
        c_prev=c;
        Q=Q+epsq;
    end
    if Q>=Qmax
        error(['snQTCnormOptFastApprox3: for T=' num2str(T) 'Q got too big']);
    end
    if T>=Tmax
        error('snQTCnormOptFastApprox3: T got too big');
    end
    if c-Kr/T >= T_limit_st
        break;
    end
    T=T+epst;
end
copt = snQTCnorm(sopt,qopt,topt,Kr,K0,L,mi,sigma,h,p);
% finally, compare with opt-base-stock
if copt > c_opt_st
    sopt = s_opt_st;
    qopt=qnot;
    topt=t_opt_st;
    copt=c_opt_st;
end

end


function [s_opt_st, t_opt_st, c_opt_st] = STCnormOptFastApprox(T_star_eoq,Kr,K0,L,mi,sigma,h,p,epst,qnot)
T=T_star_eoq;
Tmax = 100;
c_prev = +Inf;
c_opt_st = c_prev;
s_opt_st = NaN;
t_opt_st = NaN;
% smax = 10^30;
% smin = -10^30;
step_s = 0.01;
s0 = mi.*(L+T);
while T<Tmax
%     [sqt, val]=lsqnonlin(@(x) aTLCeq(x,qnot,T,L,mi,sigma,h,p), s0, smin, smax);
%     if abs(val) > 10^-4
%         error(['h=' num2str(h) ' p=' num2str(p) ' T=' num2str(T) ' s=' num2str(sqt) ' not solved val=' num2str(val)]);
%     end
    %sqt = findmins2(qnot,T,L,mi,sigma,h,p,step_s,s0);
    sqt = findmins3(qnot,T,L,mi,sigma,h,p,step_s,s0);
    s0 = sqt;
    c = snQTCnorm(sqt,qnot,T,Kr,K0,L,mi,sigma,h,p);
    %disp(['STCnormOptFastApprox(T=' num2str(T) '): sqt=' num2str(sqt) ' val=' num2str(c)]);
    if c<c_opt_st
        s_opt_st=sqt;
        t_opt_st=T;
        c_opt_st = c;
    end
    if c>c_prev
        break;
    end
    c_prev = c;
    T=T+epst;
end
if T>=Tmax
    error('STCnormOptFastApprox: T got too big');
end
end


function sqt = findmins2(Q,T,L,mi,sigma,h,phat,step_s,s0)
s = s0;
sqt = s;
c = snQTCnorm(s,Q,T,0,0,L,mi,sigma,h,phat);
cs = snQTCnorm(s+step_s,Q,T,0,0,L,mi,sigma,h,phat);
if c < cs 
   % keep decreasing s
   while true
       s = s-step_s;
       cs = snQTCnorm(s,Q,T,0,0,L,mi,sigma,h,phat);
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
            cs = snQTCnorm(s,Q,T,0,0,L,mi,sigma,h,phat);
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

