function [sopt Sopt Topt copt sopti Sopti flen Ti ci] = sSTCpoissonOpt(Tmax,Kr,K0,L,lamda,h,p,p2, epst,tnot, estop)
% the algorithm implements an approximate serial search on T.
% It uses the Zheng & Federgruen algorithm for computing optimal (s,S)
% for any given T.
% the arrays sopti and Sopti are returned because they are independent of
% Kr and therefore may subsequently be used to find solutions to other
% problems that are the same and only differ in Kr.
tic;
if nargin < 8
    p2 = 0;
end
if nargin < 9
    epst = 1;
end
if nargin < 10
    tnot = 10^-6;
end
if nargin < 11
    estop = 1;
end

copt = 10^30;

Tlen = (Tmax-tnot)/epst;
ci = 1:Tlen;
Ti = 1:Tlen;
sopti = 1:Tlen;
Sopti = 1:Tlen;
lbi = 1:Tlen;
s0 = 1;
i=1;
clprev=10^30;
for i=1:Tlen
   T = tnot + (i-1)*epst;
   Ti(i)=T;
   [s S c] = sSTCpoissonOpt_FixedT(T,Kr,K0,L,lamda,h,p,p2);
   ci(i)=c;
   sopti(i)=s;
   Sopti(i)=S;
   disp(['copt(' num2str(T) ')=' num2str(c)]);
   if c < copt
       copt = c;
       sopt = s;
       Sopt = S;
       Topt = T;
   end
   % figure out if we should stop the search
   if estop >=0
    sqt = findmins2(1,T,L,lamda,h,p,p2,s0);
    s0 = sqt;
    cl = snQTCpoisson(sqt,1,T,0,0,L,lamda,h,p,p2);  % used to be ...,Kr,0,...
    lbi(i)=cl;
    if estop == 1 && cl >= copt && cl>clprev
        break;
    end
    clprev=cl;
   end
end
tElapsed=toc;
disp(['runtime=' num2str(tElapsed) ' (secs)']);
% last stupid effort: see if (S-1,S,T) is also optimal
%c2 = sSTCpoisson(Sopt-1,Sopt,Topt,Kr,K0,L,lamda,h,p,p2);
%gap = 100*(c2-copt)/copt;
%if gap<0.1 % gap less than 0.1%
%    sopt=Sopt-1;
    %copt=c2;
%end
flen=i;
%hold on
%plot(Ti(1:i),ci(1:i));
%plot(Ti(1:i),lbi(1:i),'g');
%plot(Ti(1:i),sopti(1:i),'g');
%plot(Ti(1:i),Sopti(1:i),'y');
%hold off

end


function sqt = findmins2(Q,T,L,lamda,h,phat,p,s0)
tol=10^-5;
s = s0;
sqt = s;
c = snQTCpoisson(s,Q,T,0,0,L,lamda,h,phat,p);
cs = snQTCpoisson(s+1,Q,T,0,0,L,lamda,h,phat,p);
if c < cs 
   % keep decreasing s
   while true
       s = s-1;
       cs = snQTCpoisson(s,Q,T,0,0,L,lamda,h,phat,p);
       if cs >= c+tol  % used to be cs >= c
           sqt = s+1;
           break;
       end
       c = cs;
       sqt = s;
   end
else %  used to be else if c > cs
        % keep increasing s
        while true
            s = s+1;
            cs = snQTCpoisson(s,Q,T,0,0,L,lamda,h,phat,p);
            if cs >= c+tol  % used to be cs >= c
                sqt = s-1;
                break;
            end
            c = cs;
            sqt = s;
        end
end
end


