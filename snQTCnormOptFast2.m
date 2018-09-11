function [sopt Qopt Topt copt Qmax Tmax copt2] = snQTCnormOptFast2(Kr,K0,L,mi,sigma,h,p,epsq,epst,qnot,tnot)
copt = 10.0^30;   % the optimal (s,nQ,T) value; 
copt2 = 10.0^30;  % the optimal (s,nQ,T) value s.t. Q>epsilon
stop_at_first_Q = false; 
if nargin <11
    tnot=((3*sigma)/mi)^2;
end
if nargin <10
    qnot=10^-6;
end
if nargin < 9  % time dimension precision
    epst = 1;
end
if nargin < 8  % batch size dimension precision
    epsq = 1;
end
if nargin == 11 && tnot < 0
    stop_at_first_Q=true;
    tnot=((3*sigma)/mi)^2;
end
T=tnot;
Qmax=0;
fudge = 0.0001;  % safety feature to account for issues with the derivation of P0
while T>0
    % Q=1e-3;
    Q=qnot;
    Hmdtprev = 10^25;
    Hm = 10^25;
    popprev = 10^25;
    while Q>0
        if Q>2000
            error('Q got too big...');
        end
        % let s' be the argmin. of c(s,Q,T)
        s0=mi*(L+T);
        %smax=(mi+10.0*sigma)*(L+T);
        %smin = 0;
        %smax = (mi+50.0*sigma)*(L+T);
        %smin = -smax;
        smax = 10^30;
        smin = -10^30;
        [sqt val]=lsqnonlin(@(x) aTLCeq(x,Q,T,L,mi,sigma,h,p), s0, smin, smax);
        if abs(val) > 10^-3
            error(['h=' num2str(h) ' p=' num2str(p) ' Q=' num2str(Q) ' T=' num2str(T) ' s=' num2str(sqt) ' not solved val=' num2str(val)]);
        end
        [c Hmdt] = snQTCnormM2(sqt,Q,T,Kr,K0,L,mi,sigma,h,p);
        if c < copt
           sopt = sqt;
           Qopt = Q;
           Topt = T;
           copt = c;
        end
        if Hmdt < Hm
            Hm = Hmdt;
        end
        if c < copt2 && Q>qnot
            copt2 = c;
        end
        % find the Hmdt = min_{s}H(s,Q,T;Kr)
        %sqt2 = lsqnonlin(@(x) aTLCeq(x,Q,T,L,mi,sigma,h,p), s0, smin, smax);
        %Hmdt = snQTCnorm(sqt,Q,T,Kr,0,L,mi,sigma,h,p);
        % b = mi*T+3*sigma*sqrt(T);
        %if Hmdt + K0*mi/b > copt && Hmdt > Hmdtprev
        if Hmdt > copt && Hmdt > Hmdtprev
            break;
        end
        % break if Po(Q,T) has passed its min, and -dH/dQ is less than
        % dPo/dQ
        [pop hp passmin] = derivs(sqt,Q,T,K0,L,mi,sigma,h,p,popprev);
        if passmin==1 && pop > -hp + fudge
            break;
        end
        popprev = pop;
        Hmdtprev=Hmdt;
        Q=Q+epsq;
        if Qmax < Q
            Qmax=Q;
        end
        if stop_at_first_Q 
            break;
        end
    end
    % find the Hm = min_{s,Q}H(s,Q,T;Kr)
    %v0=[sqt2 Q];
    %vsol = fminsearch(@(v) aTLCeqH(v,T,Kr,L,mi,sigma,h,p), v0);
    %Hm = snQTCnorm(vsol(1), vsol(2), T, Kr,0,L,mi,sigma,h,p);
    % b = mi*T+3*sigma*sqrt(T);
    % if Hm + K0*mi/b > copt
    Hm2 = snQTCnorm2(sqt,qnot,T,0,0,L,mi,sigma,h,p);
    if Hm2 > copt  % used to be Hm > copt
        break;
    end
    T=T+epst;
end
Tmax=T;
end


function [Pop Hp passmin] = derivs(s,Q,T,K0,L,mi,sigma,h,p,popprev)
hq=10^-7;
l=mi;
D=sigma^2;
Pop = (K0/T)*(Porder(Q+hq,T,l,D)-Porder(Q-hq,T,l,D))/(2*hq);
if Pop > popprev
    passmin = 1;
else
    passmin = 0;
end
Hp = (snQTCnorm2(s,Q+hq,T,0,0,L,mi,sigma,h,p)-snQTCnorm2(s,Q-hq,T,0,0,L,mi,sigma,h,p))/(2*hq);
end


function y=Porder(Q,T,l,D)
P01 = l*T*(normcdf((Q-l*T)/sqrt(D*T)))/Q;
P02 = 1-normcdf((Q-l*T)/sqrt(D*T));
P03 = normpdf((Q-l*T)/sqrt(D*T))*sqrt(D*T)/Q;
y=P01+P02-P03;
end
