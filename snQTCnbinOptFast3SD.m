function [sopt Qopt Topt copt Qmax Tmax copt2] = snQTCnbinOptFast3SD(minq, mint, Kr,K0,L,lamda,pl,h,p,epsq,epst, estop, tmin, tmax)
copt = 10.0^30;   % the optimal (s,nQ,T) value; 
copt2 = 10.0^30;  % the optimal (s,nQ,T) value s.t. Q>epsilon
if nargin < 14
    tmax = 100;
end
if nargin < 13
    tmin = 5/lamda;
end
if nargin < 12  % if this value is one, run only for Q=0, i.e. run (R,T) policy
    estop = 0;
end
if nargin < 11  % time dimension precision
    epst = 1;
end
if nargin < 10  % batch size dimension precision
    epsq = 1;
end
T = tmin;
%T=1.e-6;  % start with continuous review cost (T~0)
%T = 5 / lamda;  % temporary Tmin
Qmax=0;
Qupper = 1000;
%options = optimset('MaxFunEvals',5000, 'TolFun', 1.0e-3, 'MaxIter', 3000);
s0 = 1;
s1 =1;
while T<=tmax+1.e-12
    Q=epsq;
    Hmdtprev = 10^25;
    while Q>0
        if Q>Qupper
            error('Q got too big');
        end
        %sqt = findmins2(Q,T,L,lamda,pl,h,p,0,s0);  % 0 is the p2 value
        sqt = findmins3(Q,T,L,lamda,pl,h,p,1,s0);  % 1 is eps_s
        if Q==epsq
            s1 = sqt;
        end
        s0 = sqt;
        c = snQTCnbin(sqt,Q,T,Kr,K0,L,lamda,pl,h,p);
        %disp(['for Q=' num2str(Q) ',T=' num2str(T) ' bestcost=' num2str(c)]);
        %if exitflag<=0
        %    error('optimization in snQTpoissonOptFast() failed');
        %end
        %sqt = fminsearch(@(x) snQTCpoisson(x,Q,T,Kr,K0,L,lamda,h,p), s0, options);
        %c = snQTCpoisson(sqt,Q,T,Kr,K0,L,lamda,h,p);
        if c < copt
           sopt = sqt;
           Qopt = Q;
           Topt = T;
           copt = c;
           %disp(['found better value ' num2str(copt) ' Q=' num2str(Q) ' T=' num2str(T)]);  % itc: HERE rm asap
        end
        if c < copt2 && Q>epsq
            copt2 = c;
        end
        % find the Hmdt = min_{s}H(s,Q,T;Kr)
        Hmdt = snQTCnbin(sqt,Q,T,Kr,0,L,lamda,pl,h,p);
        if Hmdt > copt && Hmdt > Hmdtprev && Q>minq
            break;
        end
        Hmdtprev=Hmdt;
        if estop==1
            break;  % only run for Q=1 (i.e. (R,T) policy)
        end
        Q=Q+epsq;
        if Qmax < Q
            Qmax=Q;
        end
    end
    % find the Hm = min_{s,Q}H(s,Q,T;Kr)
    %v0=[sqt Q];
    % vmin = [0 1];
    % vmax = [(L+T)*(lamda+20*sqrt(lamda)) 10^30];
    % vsol = fmincon(@(v) aTLCpoisseq(v,T,Kr,L,lamda,h,p), v0,[],[],[],[],vmin,vmax,options);
    % vsol = fminsearch(@(v) aTLCpoisseqH(v,T,Kr,L,lamda,h,p), v0, options);
    % Hm = snQTCpoisson(vsol(1), vsol(2), T, Kr,0,L,lamda,h,p);
    
    % itc20181004: Hm used to be computed via formula below, that computes
    % optimal holding+backorder costs for an (S,T) policy. However, we
    % don't know if this function is always increasing for all possible
    % distributions, and this is why we go for the EOQ lower bound computed
    % further below; notice we can't add the fixed costs to the lower bound
    % as Song does, because when fixed-costs are taken into account, (S,T)
    % is no longer optimal and thus is no longer a lower bound to the
    % (r,nQ,T).
    %Hm = snQTCnbin(s1,epsq,T,0,0,L,lamda,pl,h,p);  % used to be ...T,Kr,0,...
    % % itc20181207: commented line below because things seem way too slow
    Hm = -h*p*lamda*T*pl/(2*(h+p)*(1-pl)*log(1-pl));
    if Hm > copt && T>mint
        break;
%     else
%         disp(['Hm=' num2str(Hm) ' copt=' num2str(copt) ' T=' num2str(T)]);  % itc: HERE rm asap  
    end
    T=T+epst;
end
Tmax=T;
end

function sqt = findmins2(Q,T,L,lamda,pl,h,phat,p,s0)
s = s0;
sqt = s;
c = snQTCnbin(s,Q,T,0,0,L,lamda,pl,h,phat);
cs = snQTCnbin(s+1,Q,T,0,0,L,lamda,pl,h,phat);
if c < cs 
   % keep decreasing s
   while true
       s = s-1;
       cs = snQTCnbin(s,Q,T,0,0,L,lamda,pl,h,phat);
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
            cs = snQTCnbin(s,Q,T,0,0,L,lamda,pl,h,phat);
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

