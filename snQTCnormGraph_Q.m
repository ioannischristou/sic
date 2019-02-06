function [Qopt costs] = snQTCnormGraph_Q(s,Qmax,T,Kr,K0,L,mi,sigma,h,p, epsq)
% plot the c(s,Q,T) cost of the (s,nQ,T) policy as a function of Q
% for given s,T and compare with (R,T) policy

if nargin < 11
    epsq = 1.0;
end

% first figure out the cost of the (S,T) policy for the given T
[sqt tunused c_ST] = STCnormOpt2(Kr,K0,L,mi,sigma,h,p,T);

Qlen = Qmax/epsq;

costs=1:Qlen;
Qi=1:Qlen;
hmax_Qi=1:Qlen;  % h(Q) -> the min_{s} val of the cost function H(s,Q,T;Kr+K0)
hmin_Qi=1:Qlen;  % h(T) -> as above but with Kr only
c_STi = 1:Qlen;

Qopt = -1;
copt=10.0^30;
for Q=1:Qlen
    Qi(Q)=(Q-1)*epsq + 1.e-6;  % 1.e-xxx is substitute for Q=0.
    sqt = s;
    c = snQTCnorm2(sqt,Qi(Q),T,Kr,K0,L,mi,sigma,h,p);
    hmax_Qi(Q) = snQTCnorm2(sqt,Qi(Q),T,Kr+K0,0,L,mi,sigma,h,p);
    hmin_Qi(Q) = snQTCnorm2(sqt,Qi(Q),T,Kr,0,L,mi,sigma,h,p);
    costs(Q)=c;
    if c < copt
       copt = c;
       Qopt=Qi(Q);
    end
end
hold on
plot(Qi(2:Qlen),hmin_Qi(2:Qlen),'g-');
hold off
hold on
plot(Qi(2:Qlen), hmax_Qi(2:Qlen), 'r-');
hold off
for i=1:Qlen
    c_STi(i)=c_ST;
end
hold on
plot(Qi,c_STi,'b-');
hold off
hold on
plot(Qi(2:Qlen),costs(2:Qlen),'k.-');
hold off
end

