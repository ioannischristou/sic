% [sopt Qopt costs]=snQTCnormOptGraph_Q2(Qmax,T,Kr,K0,L,mi,sigma,h,p,espq)
% plot the min_{s}c(s,Q,T) cost of the (s,nQ,T) policy as a function of Q
% for given T, and compare with (R,T) policy
% if T<0, then plot the function c(s*,b(T),T), T=1...Qmax
function [sopt Qopt costs] = snQTCnormOptGraph_Q2(Qmax,T,Kr,K0,L,mi,sigma,h,p, epsq, drawstart)

if nargin < 11
    drawstart = 10;
end

if nargin < 10
    epsq = 1.0;
end

Targ = T;  % initial value of T passed in
if Targ<0
    T=0;
end

% first figure out the cost of the (S,T) policy for the given T
if Targ > 0
    [sqt tunused c_ST] = STCnormOpt2(Kr,K0,L,mi,sigma,h,p,T);
else
    c_ST=-1;
end

Qlen = Qmax/epsq;

costs=1:Qlen;
Qi=1:Qlen;
B1_Qi=1:Qlen;  % B1(Q) -> the min_{s} val of the cost function H(s,Q,T;Kr+K0)
B0_Qi=1:Qlen;  % h(T) -> the min_{s} val of the cost function Kr/T + K0*mi/Q + H(s,Q,T;0)
H0_Qi=1:Qlen;  % H(Q) -> B1(Q)-K0/T
BL_Qi=1:Qlen;  % BL(Q) -> min{B0(Q,T)-K0mi(b-Q)/bQ, B1(Q,T)-K0(b-miT)/bT, B1(a,T)-K0(1-P0(Q,T))/T}
c_STi = 1:Qlen;

a = mi*T - 3*sigma*sqrt(T);
b = mi*T + 3*sigma*sqrt(T);

Qopt = -1;
copt=10.0^30;
B0opt=10.0^30;
B0optarg=0;
first_time=1; B1_Qia = 10.0^30;
Qopttr=-1;
copttr=-10.0^30;
start = 1;
for Q=1:Qlen
    Qi(Q)=(Q-1)*epsq + 1.e-6;  % 1.e-xxx is substitute for Q=0.
    if Targ<0
        T=Qi(Q)/mi;
        if T<(3*sigma/mi)^2
            start = Q;
        end
    end
    s0=mi*(L+T);
    %smin=-(mi+10.0*sigma)*(L+T);
    %smax=(mi+10.0*sigma)*(L+T);
    smin = 0;
    smax = 10^6;
    sqt=lsqnonlin(@(x) aTLCeq(x,Qi(Q),T,L,mi,sigma,h,p), s0, smin, smax);
    c = snQTCnorm(sqt,Qi(Q),T,Kr,K0,L,mi,sigma,h,p);
    B1_Qi(Q) = snQTCnorm(sqt,Qi(Q),T,Kr+K0,0,L,mi,sigma,h,p);
    B0_Qi(Q) = K0*mi/Qi(Q) + snQTCnorm(sqt,Qi(Q),T,Kr,0,L,mi,sigma,h,p);
    if B0_Qi(Q)<B0opt
        B0opt=B0_Qi(Q);
        B0optarg=Qi(Q);
    end
    if first_time==1 && Qi(Q)>=a && Targ>0
        B1_Qia = B1_Qi(Q);
        first_time=0;
    end
    H0_Qi(Q) = B1_Qi(Q)-K0/T;
    costs(Q)=c;
    if c < copt
       copt = c;
       sopt = sqt;
       Qopt=Qi(Q);
    end
    if a<=Qi(Q) && Qi(Q)<=b && Targ>0
        bb = B0_Qi(Q) - K0*mi*(b-Qi(Q))/(Qi(Q)*b);
        bg = B1_Qi(Q) - K0*(b-mi*T)/(b*T);
        be = B1_Qia + snQTCnorm(0,Qi(Q),T,0,K0,L,mi,sigma,0,0) - K0/T;
        tmp=max(bb,bg);
        BL_Qi(Q)=max(tmp,be);
        if c>copttr
            copttr=c;
            Qopttr=Qi(Q);
        end
    else
        BL_Qi(Q)=costs(Q);
    end
end
drawstart = max(start,drawstart);
hold on
plot(Qi(drawstart:Qlen),B0_Qi(drawstart:Qlen),'y-');
hold off
hold on
plot(Qi(start:Qlen), B1_Qi(start:Qlen), 'r-');
hold off
if Targ > 0
    for i=1:Qlen
        c_STi(i)=c_ST;
    end
    hold on
    plot(Qi,c_STi,'b-');
    hold off
    hold on
    plot(mi*T-3*sqrt(T)*sigma,0,'*');
    plot(mi*T+3*sqrt(T)*sigma,0,'*');
    plot(B0optarg,0,'y*');
    % plot(Qopttr,0,'k*');
    plot(Qopt,1,'ko');
    hold off
    hold on
    plot(Qi(start:Qlen),BL_Qi(start:Qlen),'g-');
    hold off
end
hold on
plot(Qi(start:Qlen),costs(start:Qlen),'k-');
hold off
end

