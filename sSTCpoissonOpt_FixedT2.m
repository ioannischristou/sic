function [sopt Sopt copt] = sSTCpoissonOpt_FixedT2(T,Kr,K0,L,lamda,h,p,p2)
% assuming the lower bound is increasing, use the Christou Algorithm 
% to find the optimal parameters (s,S) for given T
if nargin < 8
    p2 = 0;
end
Dmax = 1000*(L+T)*lamda;  % set a bound to stop the search if it goes hey-wire

copt = 10^30;
S0=1;
S2=1;
d = 1;
while true
    [s S c] = sSTCpoissonOpt_FixedDT(d,T,Kr,K0,L,lamda,h,p,p2, S0);
    S0 = S;
    if c < copt
        sopt = s;
        Sopt = S;
        copt = c;
    end
    [sl Sl cl] = sSTCpoissonOpt_FixedDT(d,T,Kr,0,L,lamda,h,p,p2, S2);
    S2 = Sl;
    if cl > copt
        break;
    end
    if d > Dmax
        error('search gone too far...');
    end
    d = d+1;
    dif = copt - cl;
    disp(['T=' num2str(T) ' opt.diff=' num2str(dif)]);
end

end
