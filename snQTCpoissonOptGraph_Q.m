function [sopt Qopt costs] = snQTCpoissonOptGraph_Q(Qmax,T,Kr,K0,L,lamda,h,p, epsq)
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

mi = ladma;
sigma = sqrt(lamda);

Qopt = -1;
copt=10.0^30;
for Q=1:Qlen
    Qi(Q)=(Q-1)*epsq + 1; 
    s0=mi*(L+T);
    smin=0.0;
    smax=(mi+10.0*sigma)*(L+T);
    %sqt=lsqnonlin(@(x) aTLCeq(x,Qi(Q),T,L,mi,sigma,h,p), s0, smin, smax);
    %c = snQTCnorm(sqt,Qi(Q),T,Kr,K0,L,mi,sigma,h,p);
    [sqt c exitflag] = fmincon(@(x) snQTCpoisson(x,Q,T,Kr,K0,L,lamda,h,p), s0, [],[],[],[], smin, smax);
    if exitflag<=0
        error('optimization in snQTpoissonOptFast() failed');
    end    
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

