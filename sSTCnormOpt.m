function [sopt Sopt Topt copt si Si flen Ti] = sSTCnormOpt(Tmax,Kr,K0,L,mi,sigma,h,p, p2, epsd, epst, tnot)
% the algorithm implements an approximate serial search on T.
% It uses the Zheng & Federgruen algorithm for computing optimal (s,S)
% for any given T.
if nargin < 9
    p2 = 0;
end
if nargin < 10
    epsd = 1;
end
if nargin < 11
    epst = 1;
end
if nargin < 12
    tnot = (3*sigma/mi)^2;
end

copt = 10^30;

Tlen = (Tmax-tnot)/epst;
flen = Tlen;
si = 1:Tlen;
Si = 1:Tlen;
ci = 1:Tlen;
Ti = 1:Tlen;
s0 = mi*(L+tnot); 
i=1;
for i=1:Tlen
   T = tnot + (i-1)*epst;
   Ti(i)=T;
   [s S c] = sSTCnormOpt_FixedT(T,Kr,K0,L,mi,sigma,h,p, p2, epsd);
   ci(i)=c;
   si(i)=s;
   Si(i)=S;
   disp(['copt(' num2str(T) ')=' num2str(c)]);
   if c < copt
       copt = c;
       sopt = s;
       Sopt = S;
       Topt = T;
   end
   % figure out if we should stop the search
   qnot = 10^-3;
   %smax=(mi+10.0*sigma)*(L+T);
   %smin = 0;
   %smax = (mi+50.0*sigma)*(L+T);
   %smin = -smax;
   smax = 10^30;
   smin = -10^30;
   [sqt val]=lsqnonlin(@(x) aTLCeq(x,qnot,T,L,mi,sigma,h,p), s0, smin, smax);
   if abs(val) > 10^-3
       error(['h=' num2str(h) ' p=' num2str(p) ' Q=' num2str(qnot) ' T=' num2str(T) ' s=' num2str(sqt) ' not solved val=' num2str(val)]);
   end
   %sqt = findmins3(qnot,T,L,mi,sigma,h,p,p2,s0,epsd);
   s0 = sqt;
   cl = snQTCnorm2(sqt,qnot,T,0,0,L,mi,sigma,h,p);  % used to be ...,T,Kr,0,...
   if cl >= copt
      flen = i;
      break;
   end
end

%hold on
%plot(Ti(1:i),ci(1:i));
%plot(Ti, ci, 'k');
%plot(Ti, si, 'g');
%plot(Ti, Si, 'b');
%hold off

end


function sqt = findmins3(Q,T,L,mi,sigma,h,phat,p, s0, epsd)
disp(['findmins3 Q=' num2str(Q) ' T=' num2str(T)]);
if nargin < 9
    s0 = round((mi*(L+T)-Q/2.0));
end
step = epsd;  % always start with step-size epsd
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
        disp(['step=' num2str(step)]);
        cnt=niterbnd;
    end
    [dir news] = detdir(s,Q,T,L,mi,sigma,h,phat,p,epsd);
    if dir==0
        sqt = news;
        break;
    end
    if dir>0
        if prevdir<0 && step>epsd
            step = step/mul;
            cnt=niterbnd;
            disp(['step=' num2str(step)]);
        end
        s = s+step;  % itc: HERE this should probably go before if stmt
    else % dir<0
        s = s-step;
        if prevdir>0 && step>epsd
            step = step/mul;
            cnt=niterbnd;
            disp(['step=' num2str(step)]);
        end
    end
    disp(['s=' num2str(s) ' dir=' num2str(dir) ' prevdir=' num2str(prevdir)]);    
    prevdir=dir;
end
end


function [dir news] = detdir(s,Q,T,L,mi,sigma,h,phat,p,epsd)
c = snQTCnorm(s,Q,T,0,0,L,mi,sigma,h,phat,p);
cup=snQTCnorm(s+epsd,Q,T,0,0,L,mi,sigma,h,phat,p);
if c > cup
    dir=1;
    news=s;
    return;
else if c==cup
        s2=s;
        while snQTCnorm(s2,Q,T,0,0,L,mi,sigma,h,phat,p)==c
            s2=s2+epsd;
        end
        cnew = snQTCnorm(s2,Q,T,0,0,L,mi,sigma,h,phat,p);
        if cnew < c
            dir=1;
            news=s2;
            return;
        else % cnew > c
            s2 = s;
            while snQTCnorm(s2,Q,T,0,0,L,mi,sigma,h,phat,p)==c
                s2 = s2-epsd;
            end
            cnew = snQTCnorm(s2,Q,T,0,0,L,mi,sigma,h,phat,p);
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
cdown = snQTCnorm(s-epsd,Q,T,0,0,L,mi,sigma,h,phat,p);
if c > cdown
    dir=-1;
    news=s;
    return;
else if c==cdown
        s2=s;
        while snQTCnorm(s2,Q,T,0,0,L,mi,sigma,h,phat,p)==c
            s2=s2-epsd;
        end
        cnew = snQTCnorm(s2,Q,T,0,0,L,mi,sigma,h,phat,p);
        if cnew < c
            dir=-1;
            news=s2;
            return;
        else % cnew > c
            s2 = s;
            while snQTCnorm(s2,Q,T,0,0,L,mi,sigma,h,phat,p)==c
                s2 = s2+epsd;
            end
            cnew = snQTCnorm(s2,Q,T,0,0,L,mi,sigma,h,phat,p);
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
