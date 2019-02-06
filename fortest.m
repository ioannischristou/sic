% Q=1+1.e-9;
% for i=1:Q
%     disp(i)
% end
% 

Q=1;
T=11.7;
Kr=5;
K0=5;
L=1;
lamda=0.5;
pl=0.1;
h=1;
p=10;

% len=30;
% smin=-5;
% si=1:len;
% ci=1:len;
% for i=1:len
%     si(i)=smin+i-1;
%     ci(i)=snQTCnbin(si(i),Q,T,Kr,K0,L,lamda,pl,h,p);
% end
% hold on
% plot(si,ci);
% hold off

sqt=lamda*(L+T);
num=50;
Qi=1:num;
ci=1:num;
hi=1:num;
for Q=1:num
        s0=sqt;
        sqt = findmins3(Q,T,L,lamda,pl,h,p,1,s0);
        [c hm] = snQTCnbin2(sqt,Q,T,Kr,K0,L,lamda,pl,h,p);
        ci(Q)=c;
        hi(Q)=hm+Kr/T;
end
hold on
plot(Qi,ci);
hold off
hold on
plot(Qi,hi);
hold off



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

    