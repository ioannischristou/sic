function [sopt qopt copt] = snQTCnbinOptFast_FixedT(T,Kr,K0,L,lamda,pl,h,p,epsq,qnot,estop)
if nargin < 9
    epsq = 1;
end
if nargin < 10
    qnot = epsq;
end
if nargin < 11
    estop = 0;
end

Q=qnot;
Qupper=1000;  % was 150
copt = 10^25;

r = -lamda/log(1-pl);
mean = pl*r/(1-pl);
s0 = round(mean*(L+T)-Q/2.0);

while Q>0
    if Q>Qupper
        error('Q got too big');
    end
    sqt = findmins3(Q,T,L,lamda,pl,h,p,0,s0);  % 0 is the p2 value
    s0 = sqt;
    [c hm] = snQTCnbin2(sqt,Q,T,Kr,K0,L,lamda,pl,h,p);
    %disp(['for Q=' num2str(Q) ',T=' num2str(T) ' bestcost=' num2str(c)]);
    if c < copt
       sopt = sqt;
       qopt = Q;
       copt = c;
       %disp(['found better value ' num2str(copt) ' Q=' num2str(qopt) ' T=' num2str(T)]);  % itc: HERE rm asap
    end
    % find the Hmdt = min_{s}H(s,Q,T;Kr)
    Hmdt = Kr/T+hm;
    if Hmdt > copt
        break;
    end
    if estop==1
        break;  % only run for Q=1 (i.e. (R,T) policy)
    end
    Q=Q+epsq;
end

end


function sqt = findmins3(Q,T,L,lamda,pl,h,phat,p,s0)
%disp(['findmins3 Q=' num2str(Q) ' T=' num2str(T)]);
if nargin < 9
    s0 = round((lamda*(L+T)-Q/2.0));
end
step = 1;  % always start with step-size 1
%itc20180910: niterbnd used to be 5
niterbnd = 3;
%niterbnd = 5;
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
    [dir news] = detdir(s,Q,T,L,lamda,pl,h,phat,p);
    if news ~= s
        disp(['s=' num2str(s) ' news=' num2str(news)]);  % itc: HERE rm asap
    end
    if dir==0
        sqt = news;
        break;
    end
    if dir>0
        if prevdir<0 && step>1
            step = step/mul;
            cnt=niterbnd;
            %disp(['step=' num2str(step)]);
        end
        s = s+step;  % itc: HERE this should probably go before if stmt
    else % dir<0
        s = s-step;
        if prevdir>0 && step>1
            step = step/mul;
            cnt=niterbnd;
            %disp(['step=' num2str(step)]);
        end
    end
    %disp(['s=' num2str(s) ' dir=' num2str(dir) ' prevdir=' num2str(prevdir)]);    
    prevdir=dir;
end
%disp(['findmins3() sqt=' num2str(sqt)]);
end

function [dir news] = detdir(s,Q,T,L,lamda,pl,h,phat,p)
c = snQTCnbin(s,Q,T,0,0,L,lamda,pl,h,phat);
cup=snQTCnbin(s+1,Q,T,0,0,L,lamda,pl,h,phat);
if c > cup
    dir=1;
    news=s;
    return;
else if c==cup
        s2=s;
        while snQTCnbin(s2,Q,T,0,0,L,lamda,pl,h,phat)==c
            s2=s2+1;
        end
        cnew = snQTCnbin(s2,Q,T,0,0,L,lamda,pl,h,phat);
        if cnew < c
            dir=1;
            news=s2;
            return;
        else % cnew > c
            s2 = s;
            while snQTCnbin(s2,Q,T,0,0,L,lamda,pl,h,phat)==c
                s2 = s2-1;
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
cdown = snQTCnbin(s-1,Q,T,0,0,L,lamda,pl,h,phat);
if c > cdown
    dir=-1;
    news=s;
    return;
else if c==cdown
        s2=s;
        while snQTCnbin(s2,Q,T,0,0,L,lamda,pl,h,phat)==c
            s2=s2-1;
        end
        cnew = snQTCnbin(s2,Q,T,0,0,L,lamda,pl,h,phat);
        if cnew < c
            dir=-1;
            news=s2;
            return;
        else % cnew > c
            s2 = s;
            while snQTCnbin(s2,Q,T,0,0,L,lamda,pl,h,phat)==c
                s2 = s2+1;
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

