function [sopt qopt topt copt] = snQTCnbinOptFastApprox4(Kr, K0, L, lamda, pl, h, p, p2, tmin, tmax, epsq, epst)
% use the approximation for the Porder so as to obtain two convex?
% optimization problems and report the min of the two
if nargin < 12
    epst = 0;  % no quantization in T
end

if nargin < 11
    epsq = 1;  % s and Q are already quantized by the nature of the NBIN distribution
end

if nargin < 10
    tmax = 1.5;  % was 100
end

if nargin < 9
    tmin = 5/lamda;
end

if nargin < 8
    p2 = 0;
end

T0=0.5;  % was T0=1;

%v0 = [s0, T0];  % initial estimate

options = optimset('MaxFunEvals',30000, 'TolFun', 0.002, 'MaxIter', 30000, 'TolX', 0.001);

[T1, fval1, exitflag] = fmincon(@(T) f1P(T, Kr, K0, L, lamda, pl, h, p, p2, epsq), T0, [], [], [], [], tmin, tmax, [], options);
if exitflag <= 0
    error('optimization failed');
end

T0 = 0.1;  % used to be 1.0

v0 = T0;
[v2, fval2, exitflag] = fmincon(@(v) f2P(v, Kr, K0, L, lamda, pl, h, p, p2, epsq), v0, [], [], [], [], tmin, tmax, [], options);
%[v2, fval2, exitflag] = fmincon(@(v) f2P(v, Kr, K0, L, lamda, h, p), v0, [], [], [], [], [0, 1, 0], [], [], options);
if exitflag <= 0
    error('optimization failed');
end

if epst == 0
    if fval1 < fval2
        sopt = findmins3(epsq,T1,L,lamda,pl,h,p,p2,0);
        qopt = epsq;
        topt = T1;
        copt = snQTCnbin(sopt, qopt, topt, Kr, K0, L, lamda, pl, h, p);
    else
        topt = v2;
        [sopt qopt copt] = snQTCnbinOptFast_FixedT(topt,Kr,K0,L,lamda,pl,h,p, epsq);
    end
else  % treat T quantized case (epst > 0) - Q (and s) are already quantized
    tnew1 = nearest(T1, epst);
    qnew1 = epsq;
    snew1 = findmins3(epsq,tnew1,L,lamda,pl,h,p,p2,0);
    cnew1 = snQTCnbin(snew1,qnew1,tnew1, Kr, K0, L, lamda, pl, h, p);
    tnew2 = nearest(v2, epst);
    [snew2 qnew2 cnew2] = snQTCnbinOptFast_FixedT(tnew2, Kr, K0, L, lamda, pl, h, p, epsq);
    if cnew1 < cnew2
        sopt = snew1;
        qopt = qnew1;
        topt = tnew1;
        copt = cnew1;
    else
        sopt = snew2;
        qopt = qnew2;
        topt = tnew2;
        copt = cnew2;        
    end
end

end


function y = nearest(x,eps)
if eps==0
    y=x;
    return;
end
d = floor(x/eps);
y = eps*d;
y2 = eps*(d+1);
if abs(y2-x)<abs(y-x)
    y = y2;
end
if y==0
    y = eps;
end
end


function z = f2P(v, Kr, K0, L, lamda, pl, h, p, p2, epsq)
T = v;
% compute optimal sT,QT for the B_L(s,Q,T) function
Q=epsq;
Qupper=1000;
Hmdtprev = 10^25;
copt = 10^25;

r = -lamda/log(1-pl);  % T is missing as from z computation, so we're ok
mean = pl*r/(1-pl);
s0 = round(mean*(L+T)-Q/2.0);

while Q>0
    if Q>Qupper
        error('Q got too big');
    end
    sqt = findmins3(Q,T,L,lamda,pl,h,p,0,s0);  % 0 is the p2 value
    s0 = sqt;
    Hmdt = snQTCnbin(sqt,Q,T,Kr,0,L,lamda,pl,h,p);
    c = Hmdt + K0*mean/Q;
    %disp(['for Q=' num2str(Q) ',T=' num2str(T) ' bestcost=' num2str(c)]);
    if c < copt
       sT = sqt;
       QT = Q;
       copt = c;
       disp(['approx4: found better value ' num2str(copt) ' Q=' num2str(QT) ' T=' num2str(T)]);  % itc: HERE rm asap
    end
    if Hmdt > copt && Hmdt > Hmdtprev
        break;
    end
    Hmdtprev=Hmdt;
    Q=Q+epsq;
end

z = copt; 
end


function y=f1P(T,Kr,K0,L,lamda,pl,h,p,p2,epsq)
% (nQ,r,T) model in continuous time T

Q=epsq;
s=findmins3(Q,T,L,lamda,pl,h,p,p2);

%y = snQTCnbin(s,Q,T,Kr,K0,L,lamda,pl,h,p);  
% the probability should be computed directly as above, 
% because, Q=1 does not imply Porder=1.
y = (Kr+K0)/T + snQTCnbin(s,Q,T,0,0,L,lamda,pl,h,p);
end


function sqt = findmins3(Q,T,L,lamda,pl,h,phat,p,s0)
disp(['findmins3 Q=' num2str(Q) ' T=' num2str(T)]);
if nargin < 9
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
disp(['findmins3() sqt=' num2str(sqt)]);
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

