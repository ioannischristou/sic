function [sopt qopt topt copt] = snQTCnbinOptFastApprox13(Kr, K0, L, lamda, pl, h, p, epst, tnot)
% the function computes an approximate optimum point of the (r,nQ,T) policy
% when demands are distributed according to the Negative Binomial according
% to Lagodimos.
% Notice that both s and Q are quantized, with step=1, as Nbin is a
% discrete distribution.

% missing arguments
if nargin < 9
    tnot = 5./lamda;
end
if nargin < 8
    epst = 0.01;
end

%initialize return variables
sopt=NaN;
qopt=NaN;
topt=NaN;
copt=+Inf;

% compute EOQ quantities
mtdem = -lamda.*pl/((1-pl).*log(1-pl)); 
T_star_eoq = sqrt(2*K0.*(p+h)./(mtdem.*h.*p));
Q_star_eoq = sqrt(2*K0.*mtdem.*(p+h)./(h.*p));

step_s = 1;
epsq = 1;
qnot = 1;

% compute optimal (S,T) 
start_T = max(ceil(T_star_eoq/epst)*epst,tnot);  % used to be T_star_eoq
[s_opt_st, t_opt_st, c_opt_st] = STCnbinOptFastApprox(start_T,Kr,K0,L,lamda,pl,h,p,epst,qnot);


Tmax = 100;  % max limit for T-search
Qmax = 2000; % max limit for Q-search

T=tnot;
T_limit_st = c_opt_st-Kr./t_opt_st;

% outer-loop: T-search
while T<Tmax
    Q=ceil(Q_star_eoq);
    % inner-loop: Q-search
    s0=round(mtdem*(L+T));    
    c_prev = +Inf;
    while Q<Qmax
        sqt = findmins3(Q,T,L,lamda,pl,h,p,step_s,s0);
        % sqt is the current cost minimizer given (Q,T)
        c = snQTCnbinApprox2(sqt,Q,T,Kr,K0,L,lamda,pl,h,p);
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
        error(['snQTCnbinOptFastApprox13: for T=' num2str(T) 'Q got too big']);
    end
    if T>=Tmax
        error('snQTCnbinOptFastApprox13: T got too big');
    end
    if c-Kr/T >= T_limit_st
        break;
    end
    T=T+epst;
end
copt = snQTCnbin(sopt,qopt,topt,Kr,K0,L,lamda,pl,h,p);
% finally, compare with opt-base-stock
if copt > c_opt_st
    sopt = s_opt_st;
    qopt=qnot;
    topt=t_opt_st;
    copt=c_opt_st;
end

end


function [s_opt_st, t_opt_st, c_opt_st] = STCnbinOptFastApprox(T_star_eoq,Kr,K0,L,lamda,pl,h,p,epst,qnot)
T=T_star_eoq;
Tmax = 100;
c_opt_st = +Inf;
s_opt_st = NaN;
t_opt_st = NaN;
step_s = 1;
mtdem = -lamda.*pl/((1-pl).*log(1-pl));
s0 = round(mtdem.*(L+T));
hmin = +Inf;
while T<Tmax
    sqt = findmins3(qnot,T,L,lamda,pl,h,p,step_s,s0);
    s0 = sqt;
    [c hm] = snQTCnbin2(sqt,qnot,T,Kr,K0,L,lamda,pl,h,p);
    %disp(['STCnbinOptFastApprox(T=' num2str(T) '): sqt=' num2str(sqt) ' val=' num2str(c)]);
    if c<c_opt_st
        s_opt_st=sqt;
        t_opt_st=T;
        c_opt_st = c;
    end
    if hmin > hm
        hmin = hm;
    end
    if hm>=c_opt_st
        break;
    end
    T=T+epst;
end
if T>=Tmax
    error('STCnbinOptFastApprox: T got too big');
end
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

