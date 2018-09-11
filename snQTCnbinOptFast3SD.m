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
        sqt = findmins2(Q,T,L,lamda,pl,h,p,0,s0);  % 0 is the p2 value
        if Q==epsq
            s1 = sqt;
        end
        s0 = sqt;
        c = snQTCnbin(sqt,Q,T,Kr,K0,L,lamda,pl,h,p);
        disp(['for Q=' num2str(Q) ',T=' num2str(T) ' bestcost=' num2str(c)]);
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
           disp(['found better value ' num2str(copt) ' Q=' num2str(Q) ' T=' num2str(T)]);  % itc: HERE rm asap
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
    Hm = snQTCnbin(s1,epsq,T,0,0,L,lamda,pl,h,p);  % used to be ...T,Kr,0,...
    if Hm > copt && T>mint
        break;
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

