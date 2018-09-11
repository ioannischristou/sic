function [sopt Sopt copt] = sSTCpoissonOpt_FixedTSlow(T,Kr,K0,L,lamda,h,p,p2,Dmax)
% slow grid search method on [s X S] grid.
if nargin < 8
    p2 = 0;
end
if nargin < 9
    Dmax = 100*(L+T)*lamda;
end

copt = 10^30;
S0=1;
for d=1:Dmax
    [s S c] = sSTCpoissonOpt_FixedDT(d,T,Kr,K0,L,lamda,h,p,p2, S0);
    S0 = S;
    if c < copt
        sopt = s;
        Sopt = S;
        copt = c;
    end
end

