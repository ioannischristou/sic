function [sopt Qopt costs] = snQTCpoissonOptGraph_Q2(Qmax,T,Kr,K0,L,lamda,h,p, epsq)
% plot the min_{s}c(s,Q,T) cost of the (s,nQ,T) policy as a function of Q
% for given T, and compare with (R,T) policy

if nargin < 9
    epsq = 1.0;
end

% first figure out the cost of the (S,T) policy for the given T
% [sqt tunused c_ST] = STCpoissonOpt2(Kr,K0,L,lamda,h,p,T);

Qlen = Qmax/epsq;

costs=1:Qlen;
Qi=1:Qlen;
hmax_Qi=1:Qlen;  % h(Q) -> the min_{s} val of the cost function H(s,Q,T;Kr+K0)
hmin_Qi=1:Qlen;  % h(T) -> as above but with Kr only
%c_STi = 1:Qlen;

mi = lamda;
%sigma = sqrt(lamda);

Qopt = -1;
copt=10.0^30;
for Q=1:Qlen
    Qi(Q)=(Q-1)*epsq + epsq; 
    s0=mi*(L+T);
    smin=-ceil(s0);
    % compute optimal s
    % sqt = findmins(Qi(Q),T,L,lamda,h,p,smin);
    sqt = findmins2(Qi(Q),T,L,lamda,h,p,0,smin);
    c = snQTCpoisson(sqt,Qi(Q),T,Kr,K0,L,lamda,h,p);
    hmax_Qi(Q) = snQTCpoisson(sqt,Qi(Q),T,Kr+K0,0,L,lamda,h,p);
    hmin_Qi(Q) = snQTCpoisson(sqt,Qi(Q),T,Kr,0,L,lamda,h,p);
    costs(Q)=c;
    if c < copt
       copt = c;
       sopt = sqt;
       Qopt=Qi(Q);
    end
end
hold on
plot(Qi,hmin_Qi,'g-');
hold off
hold on
plot(Qi(1:Qlen), hmax_Qi(1:Qlen), 'r-');
hold off
%for i=1:Qlen
%    c_STi(i)=c_ST;
%end
%hold on
%plot(Qi,c_STi,'b-');
%hold off
hold on
plot(Qi,costs,'k.-');
hold off
end

function [sqt c] = findmins(Q,T,L,lamda,h,p,smin)
sqt = smin;
c = snQTCpoisson(sqt,Q,T,0,0,L,lamda,h,p);
s = sqt;
while true 
    s = s+1;
    cs = snQTCpoisson(s,Q,T,0,0,L,lamda,h,p);
    if cs>=c
        break;
    end
    c = cs;
    sqt = s;
end
end

function sqt = findmins2(Q,T,L,lamda,h,phat,p,s0)
s = s0;
sqt = s;
c = snQTCpoisson(s,Q,T,0,0,L,lamda,h,phat,p);
cs = snQTCpoisson(s+1,Q,T,0,0,L,lamda,h,phat,p);
if c < cs 
   % keep decreasing s
   while true
       s = s-1;
       cs = snQTCpoisson(s,Q,T,0,0,L,lamda,h,phat,p);
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
            cs = snQTCpoisson(s,Q,T,0,0,L,lamda,h,phat,p);
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