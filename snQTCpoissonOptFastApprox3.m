function [sopt qopt topt copt sopt2 qopt2 topt2 copt2] = snQTCpoissonOptFastApprox3(Kr, K0, L, lamda, h, p, p2, tmin, tmax)
% use the approximation for the Porder so as to obtain two convex
% optimization problems and report the min of the two

if nargin < 9
    tmax = 10;  % was 100
end

if nargin < 8
    tmin = 5/lamda;
end

if nargin < 7
    p2 = 0;
end

%s0 = lamda*(L+T0);
Q0 = 1;
T0 = 1;

Qupper=150;  % was 1000

%v0 = [s0, T0];  % initial estimate

options = optimset('MaxFunEvals',30000, 'TolFun', 1.0e-5, 'MaxIter', 30000);

[T1, fval1, exitflag] = fmincon(@(T) f1P(T, Kr, K0, L, lamda, h, p, p2), T0, [], [], [], [], tmin, tmax, [], options);
if exitflag <= 0
    error('optimization failed');
end

v0 = [Q0, T0];
% surprisingly, it is this optimization that turns out to have numerical
% problems and for this reason, the weird lower bounds that should be
% simply [-10^6, 1, 0]
[v2, fval2, exitflag] = fmincon(@(v) f2P(v, Kr, K0, L, lamda, h, p, p2), v0, [], [], [], [], [1, tmin], [Qupper, tmax], [], options);
%[v2, fval2, exitflag] = fmincon(@(v) f2P(v, Kr, K0, L, lamda, h, p), v0, [], [], [], [], [0, 1, 0], [], [], options);
if exitflag <= 0
    error('optimization failed');
end

if fval1 < fval2
    sopt = findmins3(1,T1,L,lamda,h,p,p2,0);
    qopt = 1;
    topt = T1;
    copt = snQTCpoisson(sopt, qopt, topt, Kr, K0, L, lamda, h, p, p2);
    %sopt2 = v2(1);
    qopt2 = round(v2(1));  % it's still >=1
    topt2 = v2(2);
    sopt2 = findmins3(qopt2,topt2,L,lamda,h,p,0);
    copt2 = snQTCpoisson(sopt2, qopt2, topt2, Kr, K0, L, lamda, h, p, p2);
else
    %sopt = v2(1);
    qopt = round(v2(1));  % it's still >=1
    topt = v2(2);
    sopt = findmins3(qopt,topt,L,lamda,h,p,p2);
    copt = snQTCpoisson(sopt, qopt, topt, Kr, K0, L, lamda, h, p, p2);
    sopt2 = findmins3(1,T1,L,lamda,h,p,p2);
    qopt2 = 1;
    topt2 = T1;
    copt2 = snQTCpoisson(sopt2, qopt2, topt2, Kr, K0, L, lamda, h, p, p2);
end
end


function z = f2P(v, Kr, K0, L, lamda, h, p, p2)
Q=v(1);  %should be Q = round(v(1));
T = v(2);
s = findmins3(Q,T,L,lamda,h,p,p2);
z1 = snQTCpoisson(s,Q,T,0,0,L,lamda,h,p,p2);
z = Kr/T + K0*lamda/Q + z1;
end


function y=f1P(T,Kr,K0,L,lamda,h,p,p2)
% (nQ,r,T) model in continuous time T
% uses the Hadley-Whittin formulas 5-20, 5-22...24
% 5-27...29, and 5-33.

Q=1;
s=findmins3(Q,T,L,lamda,h,p,p2);

%y = snQTCpoisson(s,Q,T,Kr,K0,L,lamda,h,p);  % let the probability be computed directly
y = (Kr+K0)/T + snQTCpoisson(s,Q,T,0,0,L,lamda,h,p,p2);
end


function sqt = findmins3(Q,T,L,lamda,h,phat,p,s0)
disp(['findmins3 Q=' num2str(Q) ' T=' num2str(T)]);
if nargin < 8
    s0 = round((lamda*(L+T)-Q/2.0));
end
step = 1;  % always start with step-size 1
niterbnd = 5;
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
    [dir news] = detdir(s,Q,T,L,lamda,h,phat,p);
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
end

function [dir news] = detdir(s,Q,T,L,lamda,h,phat,p)
c = snQTCpoisson(s,Q,T,0,0,L,lamda,h,phat,p);
cup=snQTCpoisson(s+1,Q,T,0,0,L,lamda,h,phat,p);
if c > cup
    dir=1;
    news=s;
    return;
else if c==cup
        s2=s;
        while snQTCpoisson(s2,Q,T,0,0,L,lamda,h,phat,p)==c
            s2=s2+1;
        end
        cnew = snQTCpoisson(s2,Q,T,0,0,L,lamda,h,phat,p);
        if cnew < c
            dir=1;
            news=s2;
            return;
        else % cnew > c
            s2 = s;
            while snQTCpoisson(s2,Q,T,0,0,L,lamda,h,phat,p)==c
                s2 = s2-1;
            end
            cnew = snQTCpoisson(s2,Q,T,0,0,L,lamda,h,phat,p);
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
cdown = snQTCpoisson(s-1,Q,T,0,0,L,lamda,h,phat,p);
if c > cdown
    dir=-1;
    news=s;
    return;
else if c==cdown
        s2=s;
        while snQTCpoisson(s2,Q,T,0,0,L,lamda,h,phat,p)==c
            s2=s2-1;
        end
        cnew = snQTCpoisson(s2,Q,T,0,0,L,lamda,h,phat,p);
        if cnew < c
            dir=-1;
            news=s2;
            return;
        else % cnew > c
            s2 = s;
            while snQTCpoisson(s2,Q,T,0,0,L,lamda,h,phat,p)==c
                s2 = s2+1;
            end
            cnew = snQTCpoisson(s2,Q,T,0,0,L,lamda,h,phat,p);
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

